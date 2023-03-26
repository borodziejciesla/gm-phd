#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
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
        double weight = 0.0;
        StateSizeVector state = StateSizeVector::Zero();
        StateSizeMatrix covariance = StateSizeMatrix::Zero();
      };

      virtual Hypothesis PredictHypothesis(const Hypothesis & hypothesis) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;

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
        PredictBirths();
        PredictExistingTargets();
        PrepareUpdate();
      }

      void PredictBirths() {}

      void PredictExistingTargets(void) {
        // Prepare for prediction 
        PrepareTransitionMatrix();
        PrepareProcessNoiseMatrix();
        // Predict
        std::transform(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(predicted_hypothesis_),
          [this](const Hypothesis & hypothesis) {
            return PredictHypothesis(hypothesis);
          }
        );
      }
      
      void PrepareUpdate() {}

      void Update(const std::vector<Measurement> & measurements) {}

      void Prune(void) {}

      void ExtractObjects(void) {}

      double prev_timestamp_ = 0.0;
      std::vector<Object> objects_;
      std::vector<Hypothesis> hypothesis_;
      std::vector<Hypothesis> predicted_hypothesis_;
  };
};  //  namespace eot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
