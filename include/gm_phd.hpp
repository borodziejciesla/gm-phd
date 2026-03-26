#ifndef GM_PHD_INCLUDE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <execution>
#include <mutex>
#include <numbers>
#include <numeric>
#include <random>
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
    calibrations_["truncation_threshold"] = std::ref(truncation_threshold_);
    calibrations_["merging_threshold"] = std::ref(merging_threshold_);
  }

  GmPhd(const GmPhd&) = delete;
  GmPhd(GmPhd&&) = delete;
  GmPhd& operator=(const GmPhd&) = delete;
  GmPhd& operator=(GmPhd&&) = delete;
  virtual ~GmPhd(void) = default;

  /**
   * @brief Run the GmPhd filter
   *
   * @param timestamp of the current cycle, used to calculate time delta for prediction
   * @param measurements collected in the current cycle, used for measurement update
   */
  void Run(const double timestamp, const std::vector<Measurement>& measurements) {
    SetTimestamps(timestamp);
    // Run Filter
    Predict();
    Update(measurements);
    // Post Processing
    Prune();
    ExtractObjects();
    // Prepare for next cycle
    std::swap(working_hypotheses_, output_hypotheses_);
    output_hypotheses_.clear();
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

    std::for_each(std::execution::par, working_hypotheses_.begin(), working_hypotheses_.end(),
                  move);
  }

  /**
   * @brief Get the list of extracted objects from the current cycle
   *
   * @return const std::vector<Object>& list of tracked objects, where each object contains its
   * state and covariance. Note that the state is in the sensor frame, and the covariance is in the
   * state space. The list is updated at the end of each cycle, and can be used for output or
   * visualization.
   */
  const std::vector<Object>& GetObjects(void) const { return objects_; }

  double GetWeightsSum(void) const {
    return std::transform_reduce(std::execution::par, working_hypotheses_.begin(),
                                 working_hypotheses_.end(), 0.0, std::plus<>(),
                                 [](const PhdHypothesis& hypothesis) { return hypothesis.weight; });
  }

 private:
  void SetTimestamps(const double timestamp) {
    if (prev_timestamp_ != 0.0) {
      time_delta_ = timestamp - prev_timestamp_;
    }
    prev_timestamp_ = timestamp;
  }

  void Predict(void) {
    MotionModel::PrepareTransitionMatrix(time_delta_);
    MotionModel::PrepareObservationMatrix();

    PredictBirths();
    PredictExistingTargets();
  }

  void PredictBirths(void) {
    PhdHypothesis birth_hypothesis;
    birth_hypothesis.weight = 1.0;
    birth_hypothesis.state = StateVector::Ones();
    birth_hypothesis.state[1] = 0.0;
    birth_hypothesis.state[3] = 0.0;
    birth_hypothesis.covariance = StateMatrix::Identity();

    working_hypotheses_.push_back(birth_hypothesis);
  }

  void PredictExistingTargets(void) {
    // Prepare for prediction
    MotionModel::PrepareTransitionMatrix(time_delta_);

    // Predict
    auto predictor = [this](PhdHypothesis& hypothesis) {
      MotionModel::PredictHypothesis(hypothesis, ps_);
    };
    std::for_each(std::execution::par, working_hypotheses_.begin(), working_hypotheses_.end(),
                  predictor);
  }

  void Update(const std::vector<Measurement>& measurements) {
    UpdateExistedHypothesis();
    MakeMeasurementUpdate(measurements);
  }

  void UpdateExistedHypothesis(void) {
    auto update = [this](PhdHypothesis& hypothesis) {
      hypothesis.weight = (1.0 - pd_) * hypothesis.predicted_weight;
      hypothesis.state = hypothesis.predicted_state;
      hypothesis.covariance = hypothesis.predicted_covariance;
    };

    std::for_each(std::execution::par, working_hypotheses_.begin(), working_hypotheses_.end(),
                  update);
  }

  void MakeMeasurementUpdate(const std::vector<Measurement>& measurements) {
    for (const auto& z : measurements) {
      std::vector<PhdHypothesis> new_hypothseses;
      new_hypothseses.reserve(working_hypotheses_.size());

      auto create_new = [this, &z, &new_hypothseses](PhdHypothesis& hypothesis) {
        // Measure distance - if too far, skip update
        const auto distance_x =
            MotionModel::StateX(hypothesis) - MotionModel::MeasurementX(z.value);
        const auto distance_y =
            MotionModel::StateY(hypothesis) - MotionModel::MeasurementY(z.value);

        if (std::pow(distance_x, 2) + std::pow(distance_y, 2) > 25.0) {
          return;
        } else {
          // Update weight, state and covariance}
          PhdHypothesis updated_hypothesis;
          updated_hypothesis.weight =
              pd_ * hypothesis.predicted_weight *
              NormPdf<MotionModel::measurement_size>(z.value, hypothesis.predicted_measurement,
                                                     hypothesis.innovation_matrix);
          updated_hypothesis.state =
              hypothesis.predicted_state +
              hypothesis.kalman_gain * (z.value - hypothesis.predicted_measurement);
          updated_hypothesis.covariance = hypothesis.predicted_covariance_aposteriori;

          new_hypothseses.push_back(updated_hypothesis);
        }
      };

      std::for_each(working_hypotheses_.begin(), working_hypotheses_.end(), create_new);

      const auto weights_sum =
          std::transform_reduce(new_hypothseses.begin(), new_hypothseses.end(), 0.0, std::plus<>(),
                                [](const PhdHypothesis& hypothesis) { return hypothesis.weight; });
      for (auto& it : new_hypothseses) {
        it.weight /= weights_sum + kappa_;
      }
      // Add new hypotheses to the working list
      working_hypotheses_.insert(working_hypotheses_.end(), new_hypothseses.begin(),
                                 new_hypothseses.end());
    }
  }

  void Prune(void) {
    for (auto& l : working_hypotheses_) {
      if (l.weight > truncation_threshold_) {
        if (l.is_merged) {
          continue;
        }

        std::vector<PhdHypothesis> merged;
        merged.reserve(working_hypotheses_.size());

        merged.push_back(l);

        l.is_merged = true;
        const auto inverse_l_cov = l.covariance.inverse();

        auto is_close = [&l, &inverse_l_cov, this](const PhdHypothesis& n) {
          const auto diff = n.state - l.state;
          const auto mahalanobis_distance = diff.transpose() * inverse_l_cov * diff;
          return mahalanobis_distance < merging_threshold_;
        };

        for (auto& i : working_hypotheses_) {
          if (!i.is_merged && is_close(i)) {
            i.is_merged = true;
            merged.push_back(i);
          }
        }

        auto accumulated_hypothesis =
            std::accumulate(merged.begin(), merged.end(), PhdHypothesis{},
                            [](const PhdHypothesis& acc, const PhdHypothesis& n) {
                              PhdHypothesis merged = acc;
                              merged.weight += n.weight;
                              merged.state += n.weight * n.state;
                              return merged;
                            });
        if (accumulated_hypothesis.weight > 0.0) {
          accumulated_hypothesis.state /= accumulated_hypothesis.weight;
        }
        accumulated_hypothesis = std::ranges::fold_left(
            merged, accumulated_hypothesis, [](const PhdHypothesis& acc, const PhdHypothesis& n) {
              PhdHypothesis sum = acc;
              const auto diff = n.state - acc.state;
              sum.covariance +=
                  (n.weight / acc.weight) * (n.covariance + (diff) * (diff).transpose());
              return sum;
            });
        output_hypotheses_.push_back(accumulated_hypothesis);
      }
    }

    // auto truncated_hypotheses = working_hypotheses_ |
    //                             std::views::filter([this](const PhdHypothesis& n) {
    //                               return n.weight > truncation_threshold_;
    //                             }) |
    //                             std::ranges::to<std::vector<PhdHypothesis>>();

    // auto merge = [&truncated_hypotheses, this](PhdHypothesis& h) {
    //   // Select close hypotheses for merging based on Mahalanobis distance
    //   h.is_merged = true;  // Mark the current hypothesis as merged to avoid self-merging
    //   const auto inverse_h_cov = h.covariance.inverse();

    //   auto is_close = [&h, &inverse_h_cov, this](const PhdHypothesis& n) {
    //     if (n.is_merged) {
    //       return false;
    //     }
    //     const auto diff = n.state - h.state;
    //     const auto mahalanobis_distance = diff.transpose() * inverse_h_cov * diff;

    //     if (mahalanobis_distance < merging_threshold_) {
    //       n.is_merged = true;  // Mark the close hypothesis as merged to avoid multiple merging
    //     }
    //     return n.is_merged;  // Return true if the hypothesis is close and not yet merged, false
    //                          // otherwise
    //   };

    //   auto close_hypothesis = truncated_hypotheses | std::views::filter(is_close) |
    //                           std::ranges::to<std::vector<PhdHypothesis>>();

    //   // Calculate the merged hypothesis by weighted averaging
    //   auto accumulated_hypothesis = std::ranges::fold_left(
    //       close_hypothesis, PhdHypothesis{}, [](const PhdHypothesis& acc, const PhdHypothesis& n)
    //       {
    //         PhdHypothesis merged = acc;
    //         merged.weight += n.weight;
    //         merged.state += n.weight * n.state;
    //         return merged;
    //       });
    //   // Normalize the merged state by the total weight
    //   if (accumulated_hypothesis.weight > 0.0) {
    //     accumulated_hypothesis.state /= accumulated_hypothesis.weight;
    //   }
    //   // Calculate merged covariance
    //   accumulated_hypothesis = std::ranges::fold_left(
    //       close_hypothesis, accumulated_hypothesis,
    //       [](const PhdHypothesis& acc, const PhdHypothesis& n) {
    //         PhdHypothesis merged = acc;
    //         const auto diff = n.state - acc.state;
    //         merged.covariance += n.weight * (n.covariance + (diff) * (diff).transpose());
    //         return merged;
    //       });
    //   // Normalize the merged covariance by the total weight
    //   if (accumulated_hypothesis.weight > 0.0) {
    //     accumulated_hypothesis.covariance /= accumulated_hypothesis.weight;
    //   }
    //   // Add the merged hypothesis to the output list
    //   this->output_hypotheses_.push_back(accumulated_hypothesis);
    // };

    // std::for_each(truncated_hypotheses.begin(), truncated_hypotheses.end(), merge);
  }

  void ExtractObjects(void) {
    objects_.clear();
    for (const auto& hypothesis : output_hypotheses_) {
      if (hypothesis.weight > 0.5) {
        objects_.push_back({hypothesis.state, hypothesis.covariance});
      }
    }
  }

  double time_delta_ = 0.0;
  double prev_timestamp_ = 0.0;

  std::vector<Object>
      objects_;  // List of extracted objects from the current cycle, to be used for output
  std::vector<PhdHypothesis>
      working_hypotheses_;                        // Hypotheses being processed in the current cycle
  std::vector<PhdHypothesis> output_hypotheses_;  // Hypotheses generated from the current cycle, to
                                                  // be used in the next cycle

  std::mutex push_mutex_;

  double pd_ = 0.8;        // Probability of detection
  double ps_ = 0.8;        // Probability of survival
  double kappa_ = 1.0e-9;  // Clutter intensity

  double truncation_threshold_ = 0.1;
  double merging_threshold_ = 3.0;
};
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
