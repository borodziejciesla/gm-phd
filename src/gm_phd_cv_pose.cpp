#include "ghm_phd_cv_pose.hpp"

namespace mot {
  GmPhdCvPose::GmPhdCvPose(const GmPhdCalibrations<4u, 2u> & calibrations)
    : GmPhd<4u, 2u>(calibrations) {}

  GmPhdCvPose::Hypothesis GmPhdCvPose::PredictHypothesis(const Hypothesis & hypothesis) {
    static Hypothesis predicted_hypothesis;

    predicted_hypothesis.weight = calibrations_.pd * hypothesis.weight;
    predicted_hypothesis.state = transition_matrix_ * hypothesis.state;
    predicted_hypothesis.covariance = transition_matrix_ * hypothesis.covariance * transition_matrix_.transpose()
      + process_noise_covariance_matrix_;

    return predicted_hypothesis;
  }

  void GmPhdCvPose::PrepareTransitionMatrix(void) {
    transition_matrix_ = StateSizeMatrix::Zero();

    transition_matrix_(0u, 0u) = 1.0;
    transition_matrix_(0u, 2u) = time_delta;

    transition_matrix_(1u, 1u) = 1.0;
    transition_matrix_(1u, 3u) = time_delta;

    transition_matrix_(2u, 2u) = 1.0;

    transition_matrix_(3u, 3u) = 1.0;
  }

  void GmPhdCvPose::PrepareProcessNoiseMatrix(void) {
    process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

    for (auto index = 0u; index < calibrations_.process_noise_diagonal.size(); index++)
      process_noise_covariance_matrix_(index, index) = calibrations_.process_noise_diagonal.at(index);
  }
} // namespace mot
