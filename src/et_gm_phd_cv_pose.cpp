#include "et_gm_phd_cv_pose.hpp"

namespace mot {
  EtGmPhdCvPose::EtGmPhdCvPose(const GmPhdCalibrations<4u, 2u> & calibrations)
    : EtGmPhd<4u, 2u>(calibrations) {
    PredictBirths();
  }

  std::pair<Eigen::Vector<double, 4u>, Eigen::Matrix<double, 4u, 4u>> 
    EtGmPhdCvPose::PredictHypothesis(const Eigen::Vector<double, 4u> & state, const Eigen::Matrix<double, 4u, 4u> & covariance) {
    const Eigen::Vector<double, 4u> predicted_state = transition_matrix_ * state;
    const Eigen::Matrix<double, 4u, 4u> predicted_covariance = transition_matrix_ * covariance * transition_matrix_.transpose()
      + time_delta * process_noise_covariance_matrix_;

    return std::make_pair(predicted_state, predicted_covariance);
  }

  void EtGmPhdCvPose::PrepareTransitionMatrix(void) {
    transition_matrix_ = StateSizeMatrix::Zero();

    transition_matrix_(0u, 0u) = 1.0;
    transition_matrix_(0u, 2u) = time_delta;

    transition_matrix_(1u, 1u) = 1.0;
    transition_matrix_(1u, 3u) = time_delta;

    transition_matrix_(2u, 2u) = 1.0;

    transition_matrix_(3u, 3u) = 1.0;
  }

  void EtGmPhdCvPose::PrepareProcessNoiseMatrix(void) {
    process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

    for (auto index = 0u; index < calibrations_.process_noise_diagonal.size(); index++)
      process_noise_covariance_matrix_(index, index) = calibrations_.process_noise_diagonal.at(index);
  }
} // namespace mot
