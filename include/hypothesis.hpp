#ifndef GM_PHD_INCLUDE_HYPOTHESIS_HPP_
#define GM_PHD_INCLUDE_HYPOTHESIS_HPP_

#include "helpers.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
struct Hypothesis {
  using StateVector = Vector<state_size>;
  using StateMatrix = Matrix<state_size>;
  using MeasurementVector = Vector<measurement_size>;
  using MeasurementMatrix = Matrix<measurement_size>;
  using KalmanGain = KalmanGainMatrix<state_size, measurement_size>;

  double weight = 0.0;
  StateVector state = StateVector::Zero();
  StateMatrix covariance = StateMatrix::Zero();

  double predicted_weight = 0.0;
  StateVector predicted_state = StateVector::Zero();
  StateMatrix predicted_covariance = StateMatrix::Zero();
  MeasurementVector predicted_measurement = MeasurementVector::Zero();
  MeasurementMatrix innovation_matrix = MeasurementMatrix::Zero();
  KalmanGain kalman_gain = KalmanGain::Zero();
  StateMatrix predicted_covariance_aposteriori = StateMatrix::Zero();

  bool is_merged = false;
};
}  // namespace mot
#endif  //  GM_PHD_INCLUDE_HYPOTHESIS_HPP_
