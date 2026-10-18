#ifndef CV_GM_PHD_HPP_
#define CV_GM_PHD_HPP_

#include "birth_model_base.hpp"
#include "gm_phd.hpp"

namespace gm_phd_example {
/* Constant Velocity (CV) motion model for the GmPhd filter */
class CvMotionModel {
 public:
  static constexpr size_t state_size = 4u;
  static constexpr size_t measurement_size = 2u;

  using StateVector = mot::Vector<state_size>;
  using StateMatrix = mot::Matrix<state_size>;
  using MeasurementVector = mot::Vector<measurement_size>;
  using MeasurementMatrix = mot::Matrix<measurement_size>;
  using Observer = mot::ObservationMatrix<state_size, measurement_size>;

  using CvHypothesis = mot::Hypothesis<state_size, measurement_size>;

  static void PrepareTransitionMatrix(const double time_delta);
  static void PrepareObservationMatrix(void);
  static void PredictHypothesis(CvHypothesis& hypothesis, const double ps);

  static double& StateX(CvHypothesis& hypothesis);
  static double& StateY(CvHypothesis& hypothesis);
  static double& StateYaw(CvHypothesis& hypothesis);

  static const double& MeasurementX(const MeasurementVector& measurement);
  static const double& MeasurementY(const MeasurementVector& measurement);

  static const Observer& GetObservationMatrix(void) { return observation_matrix_; }

 private:
  inline static StateMatrix transition_matrix_ = StateMatrix::Zero();
  inline static StateMatrix process_noise_covariance_matrix_ = StateMatrix::Zero();
  inline static MeasurementMatrix observation_noise_covariance_matrix_ =
      0.36 * MeasurementMatrix::Identity();
  inline static Observer observation_matrix_ = Observer::Zero();
};

/* Birth model for the GmPhd filter using the CV motion model */
class CvBirthModel : public mot::BirthModelBase<CvMotionModel> {
 protected:
  ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement);
  ModelHypothesis PredictBirthHypothesis(void);
  size_t PredictBirthHypothesesCount(void);

  static constexpr double measurement_birth_weight_ = 0.001;
  std::vector<ModelHypothesis> birth_hypotheses_;
};

/* Type alias for the GmPhd filter using the CV motion model and CV birth model */
using GmPhdCv = mot::GmPhd<CvMotionModel, CvBirthModel>;
using Measurement = GmPhdCv::Measurement;
}  // namespace gm_phd_example

#endif  // CV_GM_PHD_HPP