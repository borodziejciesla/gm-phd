#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include <Eigen/Dense>

#include "gm_phd_calibrations.hpp"
#include "value_with_covariance.hpp"

namespace mot {
  template <size_t state_size, size_t measurement_size>
  class GmPhd {
    public:
      using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Object = ValueWithCovariance<state_size>;
      using Measurement = ValueWithCovariance<measurement_size>;

    public:
      explicit GmPhd(const GmPhdCalibrations<state_size, measurement_size> & calibrations)
        : calibrations_{calibrations} {}

      virtual ~GmPhd(void) = default;

      void Run(const double timestamp, const std::vector<Measurement> & measurements) {
        SetTimestamps(timestamp);
        // Run Filter
        Predict();
        Update(measurements);
        // Post Processing
        Prune();
        ExtractObjects();
      }

      const std::vector<Object> & GetObjects(void) const {
        return objects_;
      }

    protected:
      struct Hypothesis {
        Hypothesis(void) = default;
        Hypothesis(const Hypothesis&) = default;
        Hypothesis(Hypothesis&&) = default;
        Hypothesis & operator=(const Hypothesis&) = default;
        Hypothesis(const double w, const StateSizeVector s, const StateSizeMatrix c)
          : weight{w}
          , state{s}
          , covariance{c} {}

        bool operator==(const Hypothesis & arg) {
          return (weight == arg.weight)
            && (state == arg.state)
            && (covariance == arg.covariance);
        }

        double weight = 0.0;
        StateSizeVector state = StateSizeVector::Zero();
        StateSizeMatrix covariance = StateSizeMatrix::Zero();
      };

      struct PredictedHypothesis {
        PredictedHypothesis(void) = default;
        PredictedHypothesis(const PredictedHypothesis&) = default;
        PredictedHypothesis(PredictedHypothesis&&) = default;
        PredictedHypothesis & operator=(const PredictedHypothesis&) = default;
        PredictedHypothesis(const Hypothesis h,
          const MeasurementSizeVector pm,
          const MeasurementSizeMatrix im,
          const Eigen::Matrix<double, state_size, measurement_size> kg,
          const StateSizeMatrix uc)
          : hypothesis{h}
          , predicted_measurement{pm}
          , innovation_matrix{im}
          , kalman_gain{kg}
          , updated_covariance{uc} {}

        Hypothesis hypothesis;

        MeasurementSizeVector predicted_measurement;
        MeasurementSizeMatrix innovation_matrix;
        Eigen::Matrix<double, state_size, measurement_size> kalman_gain;
        StateSizeMatrix updated_covariance;
      };

      virtual Hypothesis PredictHypothesis(const Hypothesis & hypothesis) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;
      virtual void PredictBirths(void) = 0;

      std::vector<PredictedHypothesis> predicted_hypothesis_;

      double time_delta = 0.0;
      GmPhdCalibrations<state_size, measurement_size> calibrations_;
      StateSizeMatrix transition_matrix_ = StateSizeMatrix::Zero();
      StateSizeMatrix process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

    private:
      void SetTimestamps(const double timestamp) {
        if (prev_timestamp_ != 0.0)
          time_delta = timestamp - prev_timestamp_;
        prev_timestamp_ = timestamp;
      }

      void Predict(void) {
        predicted_hypothesis_.clear();
        PredictBirths();
        PredictExistingTargets();
      }

      void PredictExistingTargets(void) {
        // Prepare for prediction 
        PrepareTransitionMatrix();
        PrepareProcessNoiseMatrix();
        // Predict
        std::transform(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(predicted_hypothesis_),
          [this](const Hypothesis & hypothesis) {
            const auto predicted_state = PredictHypothesis(hypothesis);

            const auto predicted_measurement = calibrations_.observation_matrix * hypothesis.state;
            const auto innovation_covariance = calibrations_.measurement_covariance
              + calibrations_.observation_matrix * hypothesis.covariance * calibrations_.observation_matrix.transpose();
            const auto kalman_gain = hypothesis.covariance * calibrations_.observation_matrix.transpose()
              * innovation_covariance.inverse();
            const auto predicted_covariance = (StateSizeMatrix::Identity() - kalman_gain * calibrations_.observation_matrix)
              * hypothesis.covariance;

            return PredictedHypothesis(predicted_state, predicted_measurement, innovation_covariance, kalman_gain, predicted_covariance);
          }
        );
      }

      void Update(const std::vector<Measurement> & measurements) {
        UpdateExistedHypothesis();
        MakeMeasurementUpdate(measurements);
      }

      void UpdateExistedHypothesis(void) {
        hypothesis_.clear();
        std::transform(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          predicted_hypothesis_.begin(),
          [this](const PredictedHypothesis & hypothesis) {
            static PredictedHypothesis updated_hypothesis;

            updated_hypothesis.hypothesis.weight = (1.0 - calibrations_.pd) * hypothesis.hypothesis.weight;
            updated_hypothesis.hypothesis.state = hypothesis.hypothesis.state;
            updated_hypothesis.hypothesis.covariance = hypothesis.hypothesis.covariance;

            return updated_hypothesis;
          }
        );
      }

      void MakeMeasurementUpdate(const std::vector<Measurement> & measurements) {
        for (const auto & measurement : measurements) {
          std::vector<Hypothesis> new_hypothesis;
          for (const auto & predicted_hypothesis : predicted_hypothesis_) {
            const auto weight = calibrations_.pd * predicted_hypothesis.hypothesis.weight * NormPdf(measurement.value, predicted_hypothesis.predicted_measurement, predicted_hypothesis.innovation_matrix);
            const auto state = predicted_hypothesis.hypothesis.state + predicted_hypothesis.kalman_gain * (measurement.value - predicted_hypothesis.predicted_measurement);
            const auto covariance = predicted_hypothesis.hypothesis.covariance;

            new_hypothesis.push_back(Hypothesis(weight, state, covariance));
          }
          // Correct weights
          const auto weights_sum = std::accumulate(new_hypothesis.begin(), new_hypothesis.end(),
            0.0,
            [](double sum, const Hypothesis & curr) {
              return sum + curr.weight;
            }
          );
          // Normalize weight
          for (auto & hypothesis : new_hypothesis)
            hypothesis.weight /= (calibrations_.kappa + weights_sum);
          // Add new hypothesis to vector
          hypothesis_.insert(hypothesis_.end(), new_hypothesis.begin(), new_hypothesis.end());
        }
      }

      void Prune(void) {
        std::vector<Hypothesis> pruned_and_merged_hypothesis;
        // Select elements with weigths over turncation threshold
        auto high_weight_filter = std::views::filter([this](const Hypothesis & hypothesis){ return (hypothesis.weight >= calibrations_.truncation_threshold); });
        auto I = hypothesis_ | high_weight_filter;
        while (~std::ranges::empty(I)) {
          // Select maximum weight element
          const auto maximum_weight_hypothesis = std::max_element(I.begin(), I.end(),
            [](const Hypothesis & a, const Hypothesis & b) {
              return a.weight > b.weight;
            }
          );
          // Select hypothesis in merging threshold
          auto close_to_maximum_filter = std::views::filter(
            [maximum_weight_hypothesis](const Hypothesis & hypothesis) {
              const auto diff = hypothesis.state - maximum_weight_hypothesis->state;
              const auto distance_matrix = diff.transpose() * hypothesis.covariance.inverse() * diff;
              return distance_matrix(0);
            }
          );
          auto L = hypothesis_ | high_weight_filter | close_to_maximum_filter;
          // Calculate new merged element
          const auto merged_weight = std::accumulate(L.begin(), L.end(),
            0.0,
            [](double sum, const Hypothesis & hypothesis) {
              return sum + hypothesis.weight;
            }
          );
          
          StateSizeVector merged_state = StateSizeVector::Zero();
          for (const auto l : L)
            merged_state += (l.weight * l.state) / merged_weight;

          StateSizeMatrix merged_covariance = StateSizeMatrix::Zero();
          for (const auto l : L) {
            const auto diff = merged_state - l.state;
            merged_covariance += (l.covariance + diff * diff.transpose()) / merged_weight;
          }

          pruned_and_merged_hypothesis.push_back(Hypothesis(merged_weight, merged_state, merged_covariance));
          // Remove L from I
          for (const auto l : L)
            std::erase(hypothesis_, l);
          // Set final hypothesis
          hypothesis_ = pruned_and_merged_hypothesis;
        }
      }

      void ExtractObjects(void) {
        objects_.clear();
        for (const auto & hypothesis : hypothesis_) {
          if (hypothesis.weight > 0.5) {
            Object object;
            object.value = hypothesis.state;
            object.covariance = hypothesis.covariance;
            objects_.push_back(object);
          }
        }
      }

      static double NormPdf(const MeasurementSizeVector & z, const MeasurementSizeVector & nu, const MeasurementSizeMatrix & cov) {
        const auto diff = z - nu;
        const auto c = 1.0 / (std::sqrt(std::pow(std::numbers::pi, measurement_size) * cov.determinant()));
        const auto e = std::exp(-0.5 * diff.transpose() * cov.inverse() * diff);
        return c * e;
      }

      double prev_timestamp_ = 0.0;
      std::vector<Object> objects_;
      std::vector<Hypothesis> hypothesis_;
  };
};  //  namespace eot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
