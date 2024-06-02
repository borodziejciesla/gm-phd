#ifndef INCLUDE_HYPOTHESIS_HPP_
#define INCLUDE_HYPOTHESIS_HPP_

#include <functional>

#include "aliases.hpp"
#include "common.hpp"
#include "value_with_covariance.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class Hypothesis {
 public:
  using Measurement = ValueWithCovariance<measurement_size>;
  using PredictionFunction = std::function<void(Hypothesis<state_size, measurement_size>&, float)>;

  inline static float pd_ = 0.0f;
  inline static float ps_ = 0.0f;
  inline static Matrix<measurement_size, state_size> observation_matrix_ =
      Matrix<measurement_size, state_size>::Zero();

  Hypothesis(void) = default;
  Hypothesis(const Hypothesis<state_size, measurement_size>&) = default;
  Hypothesis(Hypothesis<state_size, measurement_size>&&) = default;
  Hypothesis<state_size, measurement_size>& operator=(
      const Hypothesis<state_size, measurement_size>&) = default;
  Hypothesis<state_size, measurement_size>& operator=(Hypothesis<state_size, measurement_size>&&) =
      default;
  Hypothesis(const float weight, const Vector<state_size> state,
             const SquareMatrix<state_size> covariance)
      : weight_{weight}, state_{state}, covariance_{covariance} {}

  bool operator==(const Hypothesis<state_size, measurement_size>& src) {
    return (weight_ == src.weight_) && (state_ == src.state_) && (covariance_ == src.covariance_);
  }

  Hypothesis<state_size, measurement_size>& operator+(
      const Hypothesis<state_size, measurement_size>& other) {
    return *this;
  }

  Hypothesis<state_size, measurement_size> operator+(const Measurement& measurement) {
    const auto weight =
        pd_ * weight_ * NormPdf(measurement.value, predicted_measurement_, predicted_covariance_);
    const auto state =
        predicted_state_ + kalman_gain_ * (measurement.value - predicted_measurement_);
    const auto covariance = predicted_covariance_;
  }

  void PredictState(const PredictionFunction& prediction_function, const float time_delta) {
    weight_ *= (1.0f - pd_);
    prediction_function(*this, time_delta);
  }

  void UpdateExisting(void) { /**/
  }

  void MoveSensor(const SensorPoseVector& sensor_pose_delta) {
    const auto dx = state_(0) - sensor_pose_delta(0);
    const auto dy = state_(1) - sensor_pose_delta(1);
    const auto cos_dyaw = std::cos(sensor_pose_delta(2));
    const auto sin_dyaw = std::sin(sensor_pose_delta(2));
    state_(0) = cos_dyaw * dx - sin_dyaw * dy;
    state_(1) = sin_dyaw * dx + cos_dyaw * dy;
  }

  float& Weight(void) { return weight_; }

  Vector<state_size>& State(void) { return state_; }

  SquareMatrix<state_size>& Covariance(void) { return covariance_; }

  Vector<measurement_size>& PredictedMeasurement(void) { return predicted_measurement_; }

 private:
  float weight_ = 0.0f;
  Vector<state_size> state_ = Vector<state_size>::Zero();
  Vector<state_size> predicted_state_ = Vector<state_size>::Zero();
  Vector<measurement_size> predicted_measurement_ = Vector<measurement_size>::Zero();
  SquareMatrix<state_size> covariance_ = SquareMatrix<state_size>::Zero();
  SquareMatrix<measurement_size> predicted_covariance_ = SquareMatrix<measurement_size>::Zero();
  Matrix<state_size, measurement_size> kalman_gain_;
};
}  //  namespace mot

#endif  //  INCLUDE_HYPOTHESIS_HPP_
