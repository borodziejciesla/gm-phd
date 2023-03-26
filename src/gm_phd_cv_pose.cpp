#include "ghm_phd_cv_pose.hpp"

#include <random>

namespace mot {
  std::random_device r;
  std::default_random_engine e(r());

  std::uniform_real_distribution<double> pose_dist(-150, 150);
  std::uniform_real_distribution<double> velocity_dist(-10, 10);

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

  void GmPhdCvPose::PredictBirths(void) {
    for (auto index = 0; index < 100u; index++) {
      PredictedHypothesis birth_hypothesis;

      birth_hypothesis.hypothesis.weight = 1.0;

      birth_hypothesis.hypothesis.state(0u) = pose_dist(e);
      birth_hypothesis.hypothesis.state(1u) = pose_dist(e);
      birth_hypothesis.hypothesis.state(2u) = velocity_dist(e);
      birth_hypothesis.hypothesis.state(3u) = velocity_dist(e);

      birth_hypothesis.hypothesis.covariance = 5.0 * StateSizeMatrix::Identity();

      predicted_hypothesis_.push_back(birth_hypothesis);
    }
  }
} // namespace mot
