#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include <Eigen/Dense>

#include "base_gm_phd.hpp"
#include "gm_phd_calibrations.hpp"
#include "value_with_covariance.hpp"

namespace mot {
  /* Hypothesis structure */
  template <size_t state_size>
  struct GmPhdHypothesis {
    GmPhdHypothesis(void) = default;
    GmPhdHypothesis(const GmPhdHypothesis&) = default;
    //GmPhdHypothesis(GmPhdHypothesis&&) = default;
    GmPhdHypothesis & operator=(const GmPhdHypothesis&) = default;
    GmPhdHypothesis(const double w, const Eigen::Vector<double, state_size> s, const Eigen::Matrix<double, state_size, state_size> c)
      : weight{w}
      , state{s}
      , covariance{c} {}

    bool operator==(const GmPhdHypothesis & arg) {
      return (weight == arg.weight)
        && (state == arg.state)
        && (covariance == arg.covariance);
    }

    double weight = 0.0;
    Eigen::Vector<double, state_size> state = Eigen::Vector<double, state_size>::Zero();
    Eigen::Matrix<double, state_size, state_size> covariance = Eigen::Matrix<double, state_size, state_size>::Zero();
  };

  /* Predicted Hypothesis structure */
  template <size_t state_size, size_t measurement_size>
  struct GmPhdPredictedHypothesis {
    GmPhdPredictedHypothesis(void) = default;
    GmPhdPredictedHypothesis(const GmPhdPredictedHypothesis&) = default;
    //GmPhdPredictedHypothesis(GmPhdPredictedHypothesis&&) = default;
    GmPhdPredictedHypothesis & operator=(const GmPhdPredictedHypothesis&) = default;
    GmPhdPredictedHypothesis(const GmPhdHypothesis<state_size> h,
      const Eigen::Vector<double, state_size> pm,
      const Eigen::Matrix<double, state_size, state_size> im,
      const Eigen::Matrix<double, state_size, measurement_size> kg,
      const Eigen::Matrix<double, state_size, state_size> uc)
      : hypothesis{h}
      , predicted_measurement{pm}
      , innovation_matrix{im}
      , kalman_gain{kg}
      , updated_covariance{uc} {}

    GmPhdHypothesis<state_size> hypothesis;

    Eigen::Vector<double, state_size> predicted_measurement;
    Eigen::Matrix<double, state_size, state_size> innovation_matrix;
    Eigen::Matrix<double, state_size, measurement_size> kalman_gain;
    Eigen::Matrix<double, state_size, state_size> updated_covariance;
  };

  /* GM-PHD Class */
  template <size_t state_size, size_t measurement_size, typename Hypothesis = GmPhdHypothesis<state_size>, typename PredictedHypothesis = GmPhdPredictedHypothesis<state_size, measurement_size>>
  class GmPhd : public BaseGmPhd<state_size, measurement_size, GmPhdHypothesis<state_size>, GmPhdPredictedHypothesis<state_size, measurement_size>>  {
    public:
      using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Object = ValueWithCovariance<state_size>;
      using Measurement = ValueWithCovariance<measurement_size>;

    public:
      explicit GmPhd(const GmPhdCalibrations<state_size, measurement_size> & calibrations)
        : BaseGmPhd<state_size, measurement_size, GmPhdHypothesis<state_size>, GmPhdPredictedHypothesis<state_size, measurement_size>>(calibrations) {}

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

    protected:
      virtual std::pair<StateSizeVector, StateSizeMatrix> PredictHypothesis(const StateSizeVector & state, const StateSizeMatrix & covariance) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;
      
      void PredictBirths(void) {
        constexpr auto birth_objects_number = 100u;
        for (auto index = 0; index < birth_objects_number; index++) {
          Hypothesis birth_hypothesis;

          birth_hypothesis.weight = 2.0 / static_cast<double>(birth_objects_number);

          birth_hypothesis.state(0u) = pose_dist(e);
          birth_hypothesis.state(1u) = pose_dist(e);
          birth_hypothesis.state(2u) = velocity_dist(e);
          birth_hypothesis.state(3u) = velocity_dist(e);

          birth_hypothesis.covariance = 1.0 * StateSizeMatrix::Identity();

          const auto predicted_measurement = calibrations_.observation_matrix * birth_hypothesis.state;
          const auto innovation_covariance = calibrations_.measurement_covariance + calibrations_.observation_matrix * birth_hypothesis.covariance * calibrations_.observation_matrix.transpose();
          const auto kalman_gain = birth_hypothesis.covariance * calibrations_.observation_matrix.transpose() * innovation_covariance.inverse();
          const auto predicted_covariance = (StateSizeMatrix::Identity() - kalman_gain * calibrations_.observation_matrix) * birth_hypothesis.covariance;

          predicted_hypothesis_.push_back(PredictedHypothesis(birth_hypothesis, predicted_measurement, innovation_covariance, kalman_gain, predicted_covariance));
        }
      }

      std::vector<PredictedHypothesis> predicted_hypothesis_;

      double time_delta = 0.0;
      GmPhdCalibrations<state_size, measurement_size> calibrations_;
      StateSizeMatrix transition_matrix_ = StateSizeMatrix::Zero();
      StateSizeMatrix process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

    private:
      void PredictExistingTargets(void) {
        // Prepare for prediction 
        PrepareTransitionMatrix();
        PrepareProcessNoiseMatrix();
        // Predict
        std::transform(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(predicted_hypothesis_),
          [this](const Hypothesis & hypothesis) {
            const auto [predicted_state, predicted_state_covariance] = PredictHypothesis(hypothesis.state, hypothesis.covariance);
            const auto predicted_weight = calibrations_.ps * hypothesis.weight;
            const Hypothesis predicted_hypothesis = Hypothesis(predicted_weight, predicted_state, predicted_state_covariance);

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
  };
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
