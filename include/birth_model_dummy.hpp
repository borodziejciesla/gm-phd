#ifndef MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_
#define MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_

#include "birth_model_base.hpp"
#include "cv_motion_model.hpp"

namespace mot {
class BirthModelDummy : public BirthModelBase<CvMotionModel> {
 protected:
  ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement) {
    ModelHypothesis birth_hypothesis;
    birth_hypothesis.weight = 0.25;
    birth_hypothesis.state(0) = measurement.value(0);
    birth_hypothesis.state(1) = measurement.value(1);
    birth_hypothesis.covariance = StateMatrix::Identity();

    return birth_hypothesis;
  }

  ModelHypothesis PredictBirthHypothesis(void) {
    ModelHypothesis birth_hypothesis;
    birth_hypothesis.weight = 1.0;
    birth_hypothesis.state = StateVector::Ones();
    birth_hypothesis.covariance = StateMatrix::Identity();

    return birth_hypothesis;
  }

  size_t PredictBirthHypothesesCount(void) { return 1u; }
};
}  // namespace mot

#endif  //  MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_
