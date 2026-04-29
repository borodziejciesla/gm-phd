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
/**
 * @brief Gaussian Mixture Probability Hypothesis Density (GmPhd) filter for
 * multi-object tracking. The GmPhd filter is a recursive Bayesian filter that
 * estimates the number and states of multiple targets in a given area based on
 * noisy measurements. The filter maintains a set of hypotheses, where each
 * hypothesis represents a potential target with its state, covariance and
 * weight. The filter consists of three main steps: prediction, update and
 * pruning. In the prediction step, the filter predicts the next state of the
 * system based on the motion model and the probability of survival. In the
 * update step, the filter updates the weights, states and covariances of the
 * hypotheses based on the measurements and the probability of detection. In the
 * pruning step, the filter merges close hypotheses and removes those with low
 * weights. The GmPhd filter can be used for various applications, such as radar
 * tracking, video surveillance and autonomous driving, where multiple targets
 * need to be tracked simultaneously in a noisy environment.
 *
 * @tparam MotionModel - The motion model used for predicting the state of the
 * targets. The motion model should define the state vector, state matrix,
 * measurement vector and measurement matrix, as well as the
 * PrepareTransitionMatrix, PrepareObservationMatrix and PredictHypothesis methods
 * for predicting the next state of a given hypothesis based on the motion model
 * and the probability of survival.
 * @tparam BirthModelType - The birth model used for generating new hypotheses
 * based on the current measurements.
 */
template <typename MotionModel, typename BirthModelType>
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
  /**
   * @brief Construct a new Gm Phd object
   */
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
  /**
   * @brief Set the Timestamps object
   *
   * @param timestamp
   */
  void SetTimestamps(const double timestamp) {
    if (prev_timestamp_ != 0.0) {
      time_delta_ = timestamp - prev_timestamp_;
    }
    prev_timestamp_ = timestamp;
  }

  /**
   * @brief Predict the next state of the system
   *
   * The prediction step consists of two parts: predicting the birth of new targets and predicting
   * the state of existing targets. The prediction of existing targets is done in parallel using the
   * motion model, which should be implemented by defining the PrepareTransitionMatrix,
   * PrepareObservationMatrix and PredictHypothesis methods. The PrepareTransitionMatrix and
   * PrepareObservationMatrix methods are responsible for preparing the transition and observation
   * matrices based on the current time delta, while the PredictHypothesis method is responsible for
   * predicting the next state of a given hypothesis based on the motion model and the probability
   * of survival. The predicted hypotheses are stored in the working list of hypotheses, which will
   * be used for the measurement update step.
   */
  void Predict(void) {
    MotionModel::PrepareTransitionMatrix(time_delta_);
    MotionModel::PrepareObservationMatrix();

    PredictBirths();
    PredictExistingTargets();
  }

  /**
   * @brief Predict the birth of new targets
   *
   * The birth model is responsible for generating new hypotheses based on the current measurements
   * and the predefined birth model. The predicted birth hypotheses are added to the working list of
   * hypotheses, and will be processed in the next steps of the filter. The birth model can be
   * customized by implementing the BirthModelType interface, which requires
   * defining the PredictBirthHypothesis and PredictBirthHypothesesCount methods. The former should
   * return a single birth hypothesis, while the latter should return the number of birth hypotheses
   * to be generated in each cycle.
   */
  void PredictBirths(void) {
    const auto birth_hypotheses = birth_model_.PreparedictBirthHypotheses();
    working_hypotheses_.insert(working_hypotheses_.end(), birth_hypotheses.begin(),
                               birth_hypotheses.end());
  }

  /**
   * @brief Predict the state of existing targets
   *
   * The prediction of existing targets is done in parallel using the motion
   * model, which should be implemented by defining the PrepareTransitionMatrix,
   * PrepareObservationMatrix and PredictHypothesis methods. The
   * PrepareTransitionMatrix and PrepareObservationMatrix methods are
   * responsible for preparing the transition and observation matrices based on the current time
   * delta, while the PredictHypothesis method is responsible for predicting the next state of a
   * given hypothesis based on the motion model and the probability of survival. The predicted
   * hypotheses are stored in the working list of hypotheses, which will be used for the measurement
   * update step.
   */
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

  /**
   * @brief Update the hypotheses based on the measurements
   *
   * The update step consists of two parts: updating the weights, states and covariances of existing
   * hypotheses, and creating new hypotheses based on the measurements. The updating of existing
   * hypotheses is done in parallel using the measurement model, which should be implemented by
   * defining the MakeMeasurementUpdate method. The creation of new hypotheses is done by
   * implementing the BirthModelType interface, which requires defining the PredictBirthHypothesis
   * and PredictBirthHypothesesCount methods.
   */
  void Update(const std::vector<Measurement>& measurements) {
    UpdateExistedHypothesis();
    MakeMeasurementUpdate(measurements);
  }

  /**
   * @brief Update the weights, states and covariances of existing hypotheses based on the predicted
   * values and the probability of detection. The weights are updated by multiplying the predicted
   * weight by the probability of missed detection (1 - pd), while the states and covariances are
   * updated by copying the predicted state and covariance. This step is done in parallel using the
   * standard library's parallel algorithms, which allows for efficient processing of a large
   * number of hypotheses. The updated hypotheses are stored in the working list of hypotheses,
   * which will be used for the measurement update step.
   */
  void UpdateExistedHypothesis(void) {
    auto update = [this](PhdHypothesis& hypothesis) {
      hypothesis.weight = (1.0 - pd_) * hypothesis.predicted_weight;
      hypothesis.state = hypothesis.predicted_state;
      hypothesis.covariance = hypothesis.predicted_covariance;
    };

    std::for_each(std::execution::par, working_hypotheses_.begin(), working_hypotheses_.end(),
                  update);
  }

  /**
   * @brief Update the hypotheses based on the measurements
   *
   * The measurement update step consists of updating the weights, states and covariances of
   * existing hypotheses based on the predicted values and the probability of detection. The updated
   * hypotheses are stored in the working list of hypotheses, which will be used for the measurement
   * update step.
   */
  void MakeMeasurementUpdate(const std::vector<Measurement>& measurements) {
    const auto working_hypothesis_begin = working_hypotheses_.begin();
    const auto working_hypothesis_end = working_hypotheses_.end();

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
          PhdHypothesis updated_hypothesis{hypothesis};
          updated_hypothesis.predicted_measurement =
              MotionModel::GetObservationMatrix() * hypothesis.predicted_state;
          updated_hypothesis.innovation_matrix =
              MotionModel::GetObservationMatrix() * hypothesis.predicted_covariance *
                  MotionModel::GetObservationMatrix().transpose() +
              z.covariance;
          updated_hypothesis.kalman_gain = hypothesis.predicted_covariance *
                                           MotionModel::GetObservationMatrix().transpose() *
                                           updated_hypothesis.innovation_matrix.inverse();
          updated_hypothesis.predicted_covariance_aposteriori =
              (StateMatrix::Identity() -
               updated_hypothesis.kalman_gain * MotionModel::GetObservationMatrix()) *
              hypothesis.predicted_covariance;

          updated_hypothesis.weight = pd_ * hypothesis.predicted_weight *
                                      NormPdf<MotionModel::measurement_size>(
                                          z.value, updated_hypothesis.predicted_measurement,
                                          updated_hypothesis.innovation_matrix);
          updated_hypothesis.state =
              hypothesis.predicted_state +
              hypothesis.kalman_gain * (z.value - hypothesis.predicted_measurement);
          updated_hypothesis.covariance = updated_hypothesis.predicted_covariance_aposteriori;

          new_hypothseses.push_back(updated_hypothesis);
        }
      };

      const auto first_new_hypothesis_iterator = new_hypothseses.end();
      std::for_each(working_hypothesis_begin, working_hypothesis_end, create_new);

      const auto weights_sum = std::transform_reduce(
          first_new_hypothesis_iterator, new_hypothseses.end(), 0.0, std::plus<>(),
          [](const PhdHypothesis& hypothesis) { return hypothesis.weight; });
      std::for_each(first_new_hypothesis_iterator, new_hypothseses.end(),
                    [weights_sum, this](PhdHypothesis& hypothesis) {
                      hypothesis.weight /= weights_sum + kappa_;
                    });
      // Add new hypotheses to the working list
      working_hypotheses_.insert(working_hypotheses_.end(), new_hypothseses.begin(),
                                 new_hypothseses.end());
    }
  }

  /**
   * @brief Prune the hypotheses
   *
   * The pruning step consists of merging close hypotheses and removing those with low weights.
   */
  void Prune(void) {
    std::vector<PhdHypothesis> merged;
    merged.reserve(working_hypotheses_.size());

    for (auto& l : working_hypotheses_) {
      merged.clear();

      if (l.weight > truncation_threshold_) {
        if (l.is_merged) {
          continue;
        }

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
  }

  /**
   * @brief Extract the list of objects from the current hypotheses
   *
   * The extraction step consists of selecting the hypotheses with weights above a certain threshold
   * and extracting their states and covariances as the list of tracked objects. The extracted
   * objects are stored in the objects_ member variable, which can be used for output or
   * visualization. The state of each object is in the sensor frame, and the covariance is in the
   * state space. The extraction step is done at the end of each cycle, after the prediction and
   * update steps, and before preparing for the next cycle.
   */
  void ExtractObjects(void) {
    objects_.clear();
    for (const auto& hypothesis : output_hypotheses_) {
      if (hypothesis.weight > 0.5) {
        objects_.push_back({hypothesis.state, hypothesis.covariance});
      }
    }
  }

  double time_delta_ =
      0.0;  ///< [s] Time delta between the current and previous cycle, used for prediction
  double prev_timestamp_ =
      0.0;  ///< [s] Timestamp of the previous cycle, used to calculate time delta

  std::vector<Object>
      objects_;  ///< [N/A] List of extracted objects from the current cycle, to be used for output
  std::vector<PhdHypothesis>
      working_hypotheses_;  ///< [N/A] Hypotheses being processed in the current cycle
  std::vector<PhdHypothesis> output_hypotheses_;  ///< [N/A] Hypotheses generated from the current
                                                  ///< cycle, to be used in the next cycle

  BirthModelType birth_model_{};  ///< [N/A] Birth model, responsible for generating new hypotheses
                                  ///< based on the current measurements

  std::mutex push_mutex_;  ///< [N/A] Mutex for synchronizing access to the working list of
                           ///< hypotheses during parallel processing

  double pd_ = 0.8;        ///< [-] Probability of detection
  double ps_ = 0.8;        ///< [-] Probability of survival
  double kappa_ = 1.0e-9;  ///< [-] Clutter intensity

  double truncation_threshold_ =
      0.1;  ///< [-] Threshold for pruning hypotheses based on their weights
  double merging_threshold_ =
      3.0;  ///< [-] Threshold for merging hypotheses based on their Mahalanobis distance
};
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_HPP_
