#include "cv_motion_model.hpp"

namespace mot {
void CvMotionModel::PrepareTransitionMatrix(const double time_delta) {
  transition_matrix_(0u, 0u) = 1.0;
  transition_matrix_(0u, 1u) = time_delta;
  transition_matrix_(1u, 1u) = 1.0;
  transition_matrix_(2u, 2u) = 1.0;
  transition_matrix_(2u, 3u) = time_delta;
  transition_matrix_(3u, 3u) = 1.0;
}

void CvMotionModel::PrepareObservationMatrix() {
  observation_matrix_(0u, 0u) = 1.0;
  observation_matrix_(1u, 2u) = 1.0;
}

void CvMotionModel::PredictHypothesis(CvHypothesis& hypothesis, const double ps) {
  // Predict weight
  hypothesis.predicted_weight = ps * hypothesis.weight;

  // Predict state
  hypothesis.predicted_state = transition_matrix_ * hypothesis.state;
  hypothesis.predicted_covariance =
      transition_matrix_ * hypothesis.covariance * transition_matrix_.transpose() +
      process_noise_covariance_matrix_;

  // Predict measurement
  hypothesis.predicted_measurement = observation_matrix_ * hypothesis.predicted_state;

  // Predict innovation matrix
  hypothesis.innovation_matrix =
      observation_noise_covariance_matrix_ +
      observation_matrix_ * hypothesis.predicted_covariance * observation_matrix_.transpose();

  // Predict Kalman gain
  hypothesis.kalman_gain = hypothesis.predicted_covariance * observation_matrix_.transpose() *
                           hypothesis.innovation_matrix.inverse();

  // Predict covariance aposteriori
  hypothesis.predicted_covariance_aposteriori =
      (StateMatrix::Identity() - hypothesis.kalman_gain * observation_matrix_) *
      hypothesis.predicted_covariance;
}
}  // namespace mot
