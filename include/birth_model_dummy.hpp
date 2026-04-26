#ifndef MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_
#define MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_

#include "birth_model_base.hpp"

namespace mot {
class BirthModelDummy : public BirthModelBase<4u, 2u> {
 protected:
  ModelHypothesis PredictBirthHypothesis(void) {
    ModelHypothesis birth_hypothesis;
    birth_hypothesis.weight = 0.5;
    birth_hypothesis.state = StateVector::Ones();
    birth_hypothesis.covariance = StateMatrix::Identity();

    return birth_hypothesis;
  }

  size_t PredictBirthHypothesesCount(void) { return 1u; }
};
}  // namespace mot

#endif  //  MOT_INCLUDE_BIRTH_MODEL_DUMMY_HPP_
