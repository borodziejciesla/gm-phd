#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>
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
      using SensorPoseVector = Eigen::Vector<double, 3u>;
      using SensorPoseMatrix = Eigen::Matrix<double, 3u, 3u>;

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

      void MoveSensor(const SensorPoseVector & sensor_pose_delta, const SensorPoseMatrix & sensor_pose_delta_covariance) {
        std::ignore = sensor_pose_delta;
        std::ignore = sensor_pose_delta_covariance;

        std::transform(hypothesis_.begin(), hypothesis_.end(),
          hypothesis_.begin(),
          [sensor_pose_delta,sensor_pose_delta_covariance](const Hypothesis & hypothesis) {
            const auto dx = hypothesis.state(0) - sensor_pose_delta(0);
            const auto dy = hypothesis.state(1) - sensor_pose_delta(1);
            const auto cos_dyaw = std::cos(sensor_pose_delta(2));
            const auto sin_dyaw = std::sin(sensor_pose_delta(2));
            hypothesis.state(0) = cos_dyaw * dx - sin_dyaw * dy;
            hypothesis.state(1) = sin_dyaw * dx + cos_dyaw * dy;
          }
        );
      }

      const std::vector<Object> & GetObjects(void) const {
        return objects_;
      }

      double GetWeightsSum(void) const {
        return std::accumulate(hypothesis_.begin(), hypothesis_.end(),
          0.0,
          [](double sum, const Hypothesis & hypothesis) {
            return sum + hypothesis.weight;
          }
        );
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
        if (is_initialized_)
          predicted_hypothesis_.clear();
        else
          is_initialized_ = true;
        // PredictBirths();
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
        //UpdateExistedHypothesis();
        MakeMeasurementUpdate(measurements);
      }

      void UpdateExistedHypothesis(void) {
        hypothesis_.clear();
        std::transform(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          predicted_hypothesis_.begin(),
          [this](const PredictedHypothesis & hypothesis) {
            PredictedHypothesis updated_hypothesis = hypothesis;

            updated_hypothesis.hypothesis.weight = (1.0 - calibrations_.pd) * hypothesis.hypothesis.weight;
            updated_hypothesis.hypothesis.state = hypothesis.hypothesis.state;
            updated_hypothesis.hypothesis.covariance = hypothesis.hypothesis.covariance;

            return updated_hypothesis;
          }
        );
      }

      void MakeMeasurementUpdate(const std::vector<Measurement> & measurements) {
        hypothesis_.clear();

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
            [this](double sum, const Hypothesis & curr) {
              return sum + curr.weight * (1.0 - calibrations_.pd);
            }
          );
          // Normalize weight
          for (auto & hypothesis : new_hypothesis)
            hypothesis.weight *= ((1.0 - calibrations_.pd) / (calibrations_.kappa + weights_sum));
          // Add new hypothesis to vector
          hypothesis_.insert(hypothesis_.end(), new_hypothesis.begin(), new_hypothesis.end());
        }
        // Add prediced previously
        std::transform(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          std::back_inserter(hypothesis_),
          [](const PredictedHypothesis & predicted_hypothesis) {
            return predicted_hypothesis.hypothesis;
          }
        );
      }

      void Prune(void) {
        // Select elements with weigths over turncation threshold
        std::vector<Hypothesis> pruned_hypothesis;
        std::copy_if(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(pruned_hypothesis),
          [this](const Hypothesis & hypothesis) {
            return hypothesis.weight >= calibrations_.truncation_threshold;
          }
        );
        std::vector<std::pair<Hypothesis, bool>> pruned_hypothesis_marked;
        std::transform(pruned_hypothesis.begin(), pruned_hypothesis.end(),
          std::back_inserter(pruned_hypothesis_marked),
          [](const Hypothesis & hypothesis) {
            return std::make_pair(hypothesis, false);
          }        
        );

        // Merge hypothesis
        std::vector<Hypothesis> merged_hypothesis;
        auto non_marked_hypothesis_counter = [](size_t sum, const std::pair<Hypothesis, bool> & markable_hypothesis) {
          return sum + (markable_hypothesis.second ? 0u : 1u);
        };
        auto non_merged_hypothesis_number = std::accumulate(pruned_hypothesis_marked.begin(), pruned_hypothesis_marked.end(), 0u, non_marked_hypothesis_counter);

        while (non_merged_hypothesis_number > 0u) {
          auto I = pruned_hypothesis_marked | std::views::filter([](const std::pair<Hypothesis, bool> & hypothesis_mark) { return !hypothesis_mark.second; });

          // Select maximum weight element
          const auto maximum_weight_hypothesis = *std::max_element(I.begin(), I.end(),
            [](const std::pair<Hypothesis, bool> & a, const std::pair<Hypothesis, bool> & b) {
              return a.first.weight < b.first.weight;
            }
          );

          // Select hypothesis in merging threshold
          auto L = pruned_hypothesis_marked | std::views::filter(
            [maximum_weight_hypothesis,this](const std::pair<Hypothesis, bool> & markable_hypothesis) {
              const auto diff = markable_hypothesis.first.state - maximum_weight_hypothesis.first.state;
              const auto distance_matrix = diff.transpose() * markable_hypothesis.first.covariance.inverse() * diff;
              return (distance_matrix(0) < calibrations_.merging_threshold) && !markable_hypothesis.second;
            }
          );

          // Calculate new merged element
          const auto merged_weight = std::accumulate(L.begin(), L.end(),
            0.0,
            [](double sum, const std::pair<Hypothesis, bool> & hypothesis) {
              return sum + hypothesis.first.weight;
            }
          );
          
          StateSizeVector merged_state = StateSizeVector::Zero();
          for (const auto l : L)
            merged_state += (l.first.weight * l.first.state) / merged_weight;

          StateSizeMatrix merged_covariance = StateSizeMatrix::Zero();
          for (const auto l : L) {
            const auto diff = merged_state - l.first.state;
            merged_covariance += (l.first.covariance + diff * diff.transpose()) / merged_weight;
          }

          merged_hypothesis.push_back(Hypothesis(merged_weight, merged_state, merged_covariance));
          // Remove L from I
          std::transform(L.begin(), L.end(),
           L.begin(),
            [](std::pair<Hypothesis, bool> & markable_hypothesis) {
              markable_hypothesis.second = true;
              return markable_hypothesis;
            }
          );
          //
          non_merged_hypothesis_number = std::accumulate(pruned_hypothesis_marked.begin(), pruned_hypothesis_marked.end(), 0u, non_marked_hypothesis_counter);
        }
        // Set final hypothesis
        hypothesis_ = merged_hypothesis;
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
      bool is_initialized_ = false;
      std::vector<Object> objects_;
      std::vector<Hypothesis> hypothesis_;
  };
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
