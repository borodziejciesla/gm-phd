#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <execution>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>
#include <vector>

#include "calibrated_object.hpp"
#include "helpers.hpp"
#include "hypothesis.hpp"

namespace mot {
template <typename MotionModel>
class GmPhd : public CalibratedObject {
 public:
  using StateVector = MotionModel::StateVector;
  using StateMatrix = MotionModel::StateMatrix;
  using MeasurementVector = MotionModel::MeasurementVector;
  using MeasurementMatrix = MotionModel::MeasurementMatrix;

  using Object = ValueWithCovariance<MotionModel::state_size>;
  using Measurement = ValueWithCovariance<MotionModel::measurement_size>;

  using PhdHypothesis = Hypothesis<MotionModel::state_size, MotionModel::measurement_size>;

 public:
  GmPhd() : CalibratedObject() {
    calibrations_["pd"] = std::ref(pd_);
    calibrations_["ps"] = std::ref(ps_);
    calibrations_["kappa"] = std::ref(kappa_);
  }

  GmPhd(const GmPhd&) = delete;
  GmPhd(GmPhd&&) = delete;
  GmPhd& operator=(const GmPhd&) = delete;
  GmPhd& operator=(GmPhd&&) = delete;
  virtual ~GmPhd(void) = default;

  void Run(const double timestamp, const std::vector<Measurement>& measurements) {
    SetTimestamps(timestamp);
    // Run Filter
    Predict();
    Update(measurements);
    // Post Processing
    // Prune();
    ExtractObjects();
  }

  void MoveSensor(const SensorPoseVector& sensor_pose_delta,
                  const SensorPoseMatrix& sensor_pose_delta_covariance) {
    auto move = [sensor_pose_delta, sensor_pose_delta_covariance](PhdHypothesis& hypothesis) {
      const auto dx = hypothesis.state(0) - sensor_pose_delta(0);
      const auto dy = hypothesis.state(1) - sensor_pose_delta(1);
      const auto cos_dyaw = std::cos(sensor_pose_delta(2));
      const auto sin_dyaw = std::sin(sensor_pose_delta(2));
      hypothesis.state(0) = cos_dyaw * dx - sin_dyaw * dy;
      hypothesis.state(1) = sin_dyaw * dx + cos_dyaw * dy;
    };

    std::for_each(std::execution::par, hypotheses_.begin(), hypotheses_.end(), move);
  }

  const std::vector<Object>& GetObjects(void) const { return objects_; }

  double GetWeightsSum(void) const {
    // accumulate
    auto add = [](double sum, const PhdHypothesis& hypothesis) { return sum + hypothesis.weight; };
    return std::transform_reduce(std::execution::par, hypotheses_.begin(), hypotheses_.end(), 0.0,
                                 add);
  }

 private:
  double time_delta_ = 0.0;
  double prev_timestamp_ = 0.0;

  void SetTimestamps(const double timestamp) {
    if (prev_timestamp_ != 0.0) {
      time_delta_ = timestamp - prev_timestamp_;
    }
    prev_timestamp_ = timestamp;
  }

  void Predict(void) {
    // TODO: Check if initialized
    // PredictBirths();
    PredictExistingTargets();
  }

  void PredictExistingTargets(void) {
    // Prepare for prediction
    MotionModel::PrepareTransitionMatrix(time_delta_);

    // Predict
    auto predictor = [this](PhdHypothesis& hypothesis) {
      MotionModel::PredictHypothesis(hypothesis, time_delta_);
    };
    std::for_each(std::execution::par, hypotheses_.begin(), hypotheses_.end(), predictor);
  }

  void Update(const std::vector<Measurement>& measurements) {
    // UpdateExistedHypothesis();
    // MakeMeasurementUpdate(measurements);
  }

  // void UpdateExistedHypothesis(void) {
  //   hypotheses_.clear();
  //   std::transform(predicted_hypotheses_.begin(), predicted_hypotheses_.end(),
  //                  predicted_hypotheses_.begin(), [this](const PredictedHypothesis& hypothesis) {
  //                    PredictedHypothesis updated_hypothesis = hypothesis;

  //                    updated_hypothesis.hypothesis.weight =
  //                        (1.0 - calibrations_.pd) * hypothesis.hypothesis.weight;
  //                    updated_hypothesis.hypothesis.state = hypothesis.hypothesis.state;
  //                    updated_hypothesis.hypothesis.covariance = hypothesis.hypothesis.covariance;

  //                    return updated_hypothesis;
  //                  });
  // }

  // void MakeMeasurementUpdate(const std::vector<Measurement>& measurements) {
  //   hypotheses_.clear();

  //   for (const auto& measurement : measurements) {
  //     std::vector<Hypothesis> new_hypothesis;
  //     for (const auto& predicted_hypothesis : predicted_hypotheses_) {
  //       const auto weight = calibrations_.pd * predicted_hypothesis.hypothesis.weight *
  //                           NormPdf(measurement.value,
  //                           predicted_hypothesis.predicted_measurement,
  //                                   predicted_hypothesis.innovation_matrix);
  //       const auto state = predicted_hypothesis.hypothesis.state +
  //                          predicted_hypothesis.kalman_gain *
  //                              (measurement.value - predicted_hypothesis.predicted_measurement);
  //       const auto covariance = predicted_hypothesis.hypothesis.covariance;

  //       new_hypothesis.push_back(Hypothesis(weight, state, covariance));
  //     }
  //     // Correct weights
  //     const auto weights_sum =
  //         std::accumulate(new_hypothesis.begin(), new_hypothesis.end(), 0.0,
  //                         [this](double sum, const Hypothesis& curr) {
  //                           return sum + curr.weight * (1.0 - calibrations_.pd);
  //                         });
  //     // Normalize weight
  //     for (auto& hypothesis : new_hypothesis)
  //       hypothesis.weight *= ((1.0 - calibrations_.pd) / (calibrations_.kappa + weights_sum));
  //     // Add new hypothesis to vector
  //     hypotheses_.insert(hypotheses_.end(), new_hypothesis.begin(), new_hypothesis.end());
  //   }
  //   // Add prediced previously
  //   std::transform(predicted_hypotheses_.begin(), predicted_hypotheses_.end(),
  //                  std::back_inserter(hypotheses_),
  //                  [](const PredictedHypothesis& predicted_hypothesis) {
  //                    return predicted_hypothesis.hypothesis;
  //                  });
  // }

  // void Prune(void) {
  //   // Select elements with weigths over turncation threshold
  //   std::vector<Hypothesis> pruned_hypothesis;
  //   std::copy_if(hypotheses_.begin(), hypotheses_.end(), std::back_inserter(pruned_hypothesis),
  //                [this](const Hypothesis& hypothesis) {
  //                  return hypothesis.weight >= calibrations_.truncation_threshold;
  //                });
  //   std::vector<std::pair<Hypothesis, bool>> pruned_hypotheses_marked;
  //   std::transform(pruned_hypothesis.begin(), pruned_hypothesis.end(),
  //                  std::back_inserter(pruned_hypotheses_marked),
  //                  [](const Hypothesis& hypothesis) { return std::make_pair(hypothesis, false);
  //                  });

  //   // Merge hypothesis
  //   std::vector<Hypothesis> merged_hypothesis;
  //   auto non_marked_hypotheses_counter =
  //       [](size_t sum, const std::pair<Hypothesis, bool>& markable_hypothesis) {
  //         return sum + (markable_hypothesis.second ? 0u : 1u);
  //       };
  //   auto non_merged_hypotheses_number =
  //       std::accumulate(pruned_hypotheses_marked.begin(), pruned_hypotheses_marked.end(), 0u,
  //                       non_marked_hypotheses_counter);

  //   while (non_merged_hypotheses_number > 0u) {
  //     auto I = pruned_hypotheses_marked |
  //              std::views::filter([](const std::pair<Hypothesis, bool>& hypotheses_mark) {
  //                return !hypotheses_mark.second;
  //              });

  //     // Select maximum weight element
  //     const auto maximum_weight_hypothesis = *std::max_element(
  //         I.begin(), I.end(),
  //         [](const std::pair<Hypothesis, bool>& a, const std::pair<Hypothesis, bool>& b) {
  //           return a.first.weight < b.first.weight;
  //         });

  //     // Select hypothesis in merging threshold
  //     auto L = pruned_hypotheses_marked |
  //              std::views::filter([maximum_weight_hypothesis,
  //                                  this](const std::pair<Hypothesis, bool>& markable_hypothesis)
  //                                  {
  //                const auto diff =
  //                    markable_hypothesis.first.state - maximum_weight_hypothesis.first.state;
  //                const auto distance_matrix =
  //                    diff.transpose() * markable_hypothesis.first.covariance.inverse() * diff;
  //                return (distance_matrix(0) < calibrations_.merging_threshold) &&
  //                       !markable_hypothesis.second;
  //              });

  //     // Calculate new merged element
  //     const auto merged_weight = std::accumulate(
  //         L.begin(), L.end(), 0.0, [](double sum, const std::pair<Hypothesis, bool>& hypothesis)
  //         {
  //           return sum + hypothesis.first.weight;
  //         });

  //     StateVector merged_state = StateVector::Zero();
  //     for (const auto l : L) merged_state += (l.first.weight * l.first.state) / merged_weight;

  //     StateMatrix merged_covariance = StateMatrix::Zero();
  //     for (const auto l : L) {
  //       const auto diff = merged_state - l.first.state;
  //       merged_covariance += (l.first.covariance + diff * diff.transpose()) / merged_weight;
  //     }

  //     merged_hypothesis.push_back(Hypothesis(merged_weight, merged_state, merged_covariance));
  //     // Remove L from I
  //     std::transform(L.begin(), L.end(), L.begin(),
  //                    [](std::pair<Hypothesis, bool>& markable_hypothesis) {
  //                      markable_hypothesis.second = true;
  //                      return markable_hypothesis;
  //                    });
  //     //
  //     non_merged_hypotheses_number =
  //         std::accumulate(pruned_hypotheses_marked.begin(), pruned_hypotheses_marked.end(), 0u,
  //                         non_marked_hypotheses_counter);
  //   }
  //   // Set final hypothesis
  //   hypotheses_ = merged_hypothesis;
  // }

  void ExtractObjects(void) {
    objects_.clear();
    for (const auto& hypothesis : hypotheses_) {
      if (hypothesis.weight > 0.5) {
        objects_.push_back({hypothesis.state, hypothesis.covariance});
      }
    }
  }

  bool is_initialized_ = false;
  std::vector<Object> objects_;
  std::vector<PhdHypothesis> hypotheses_;

  // calibrations
  std::array<double, MotionModel::state_size> process_noise_diagonal =
      {};  // Process noise covariance matrix

  Eigen::Matrix<double, MotionModel::measurement_size, MotionModel::state_size>
      observation_matrix;  // Observation matrix
  Eigen::Matrix<double, MotionModel::measurement_size, MotionModel::measurement_size>
      measurement_covariance;

  double pd_ = 0.8;  // Probability of detection
  double ps_ = 0.8;  // Probability of survival
  double kappa_ = 1.0e-9;

  double truncation_threshold = 0.1;
  double merging_threshold = 3.0;
};
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
