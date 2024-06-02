#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>
#include <vector>

#include "aliases.h"
#include "gm_phd_calibrations.hpp"
#include "hypothesis.hpp"
#include "objects_extractor_factory.hpp"
#include "value_with_covariance.hpp"

namespace mot {

// template <size_t state_size, size_t measurement_size>
// void PredictHypothesis(Hypothesis<state_size, measurement_size>& hypothesis, float time_delta_) {
//   //
// }

template <size_t state_size, size_t measurement_size>
class GmPhd {
 public:
  using StateHypothesis = Hypothesis<state_size, measurement_size>;
  using Object = ValueWithCovariance<state_size>;
  using Measurement = ValueWithCovariance<measurement_size>;

  using HypothesisVector = std::vector<StateHypothesis>;
  using MeasurementsVector = std::vector<Measurement>;
  using ObjectsVector = std::vector<Object>;

  explicit GmPhd(const GmPhdCalibrations<state_size, measurement_size>& calibrations)
      : calibrations_{calibrations},
        object_extractor_{
            ObjectsExtractorFactory<state_size, measurement_size>::CreateExtractor()} {
    StateHypothesis::pd_ = calibrations_.pd;
    StateHypothesis::ps_ = calibrations_.ps;
    StateHypothesis::observation_matrix_ = calibrations.observation_matrix;
  }

  virtual ~GmPhd(void) = default;

  void Run(const double timestamp, const MeasurementsVector& measurements) {
    SetTimestamps(timestamp);
    /* Run Filter */
    Predict();
    Update(measurements);
    /* Post Processing */
    Prune();
    Merge();
    ExtractObjects();
  }

  void MoveSensor(const SensorPoseVector& sensor_pose_delta) {
    for (auto& hypothesis : hypothesis_) {
      hypothesis.MoveSensor(sensor_pose_delta);
    }
  }

  const ObjectsVector& GetObjects(void) const { return objects_; }

  float GetWeightsSum(void) const {
    return std::accumulate(
        hypothesis_.begin(), hypothesis_.end(), 0.0f,
        [](float sum, const StateHypothesis& hypothesis) { return sum + hypothesis.Weight(); });
  }

 protected:
  // virtual StateHypothesis PredictHypothesis(const StateHypothesis& hypothesis) = 0;
  // virtual void PrepareTransitionMatrix(void) = 0;
  // virtual void PrepareProcessNoiseMatrix(void) = 0;
  // virtual void PredictBirths(void) = 0;

  HypothesisVector predicted_hypothesis_;

  float time_delta_ = 0.0f;
  GmPhdCalibrations<state_size, measurement_size> calibrations_;
  SquareMatrix<state_size> transition_matrix_ = SquareMatrix<state_size>::Zero();
  SquareMatrix<state_size> process_noise_covariance_matrix_ = SquareMatrix<state_size>::Zero();

 private:
  void SetTimestamps(const double timestamp) {
    if (prev_timestamp_ != 0.0) {
      time_delta_ = static_cast<float>(timestamp - prev_timestamp_);
    }
    prev_timestamp_ = timestamp;
  }

  void Predict(void) {
    if (is_initialized_)
      predicted_hypothesis_.clear();
    else
      is_initialized_ = true;

    // PredictBirths(); // TODO
    PredictExistingTargets();
  }

  void PredictExistingTargets(void) {
    // Prepare for prediction
    // PrepareTransitionMatrix(); // TODO
    // PrepareProcessNoiseMatrix(); // TODO
    // Predict
    std::transform(hypothesis_.begin(), hypothesis_.end(),
                   std::back_inserter(predicted_hypothesis_), [this](StateHypothesis& hypothesis) {
                     hypothesis.PredictState(calibrations_.predict_hypothesis, time_delta_);
                     return hypothesis;
                   });
  }

  void Update(const MeasurementsVector& measurements) {
    UpdateExistedHypothesis();
    MakeMeasurementUpdate(measurements);
  }

  void UpdateExistedHypothesis(void) {
    hypothesis_.clear();
    for (auto& hypothesis : hypothesis_) {
      hypothesis.UpdateExisting();
    }
  }

  void MakeMeasurementUpdate(const MeasurementsVector& measurements) { hypothesis_.clear(); }

  void Prune() {
    hypothesis_.clear();
    std::copy_if(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
                 std::back_inserter(hypothesis_), [this](StateHypothesis& hypothesis) {
                   return hypothesis.Weight() >= calibrations_.truncation_threshold;
                 });
  }

  void Merge(void) {
    //
  }

  void ExtractObjects(void) {
    objects_.clear();
    objects_ = object_extractor_->ExtractObjects(hypothesis_);
  }

  double prev_timestamp_ = 0.0;
  bool is_initialized_ = false;
  std::unique_ptr<ObjectExtractorInterface<state_size, measurement_size>> object_extractor_;
  ObjectsVector objects_;
  HypothesisVector hypothesis_;
};
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
