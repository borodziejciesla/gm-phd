#ifndef MOT_INCLUDE_CV_MOTION_MODEL_HPP_
#define MOT_INCLUDE_CV_MOTION_MODEL_HPP_

#include "helpers.hpp"
#include "hypothesis.hpp"

namespace mot {
class CvMotionModel {
 public:
  static constexpr size_t state_size = 4u;
  static constexpr size_t measurement_size = 2u;

  using StateVector = Vector<state_size>;
  using StateMatrix = Matrix<state_size>;
  using MeasurementVector = Vector<measurement_size>;
  using MeasurementMatrix = Matrix<measurement_size>;
  using Observer = ObservationMatrix<state_size, measurement_size>;

  using CvHypothesis = Hypothesis<state_size, measurement_size>;

  static void PrepareTransitionMatrix(const double time_delta);
  static void PrepareObservationMatrix(void);
  static void PredictHypothesis(CvHypothesis& hypothesis, const double ps);

  static double& StateX(CvHypothesis& hypothesis);
  static double& StateY(CvHypothesis& hypothesis);
  static double& StateYaw(CvHypothesis& hypothesis);

  static const double& MeasurementX(const MeasurementVector& measurement);
  static const double& MeasurementY(const MeasurementVector& measurement);

 private:
  inline static StateMatrix transition_matrix_ = StateMatrix::Zero();
  inline static StateMatrix process_noise_covariance_matrix_ = StateMatrix::Zero();
  inline static MeasurementMatrix observation_noise_covariance_matrix_ = MeasurementMatrix::Zero();
  inline static Observer observation_matrix_ = Observer::Zero();
};
}  // namespace mot

#endif  //  MOT_INCLUDE_CV_MOTION_MODEL_HPP_
