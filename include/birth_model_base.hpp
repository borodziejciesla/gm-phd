#ifndef MOT_INCLUDE_BIRTH_MODEL_HPP_
#define MOT_INCLUDE_BIRTH_MODEL_HPP_

#include <vector>

#include "helpers.hpp"
#include "hypothesis.hpp"

namespace mot {
template <size_t state_size_, size_t measurement_size_>
class BirthModelBase {
 public:
  static constexpr size_t state_size = state_size_;
  static constexpr size_t measurement_size = measurement_size_;

  using StateVector = Vector<state_size>;
  using StateMatrix = Matrix<state_size>;

  using ModelHypothesis = Hypothesis<state_size, measurement_size>;

  const std::vector<ModelHypothesis>& PreparedictBirthHypotheses(void) {
    birth_hypotheses_.clear();

    for (size_t i = 0u; i < PredictBirthHypothesesCount(); ++i) {
      birth_hypotheses_.emplace_back(std::move(PredictBirthHypothesis()));
    }

    return birth_hypotheses_;
  }

 protected:
  virtual ModelHypothesis PredictBirthHypothesis(void) = 0;
  virtual size_t PredictBirthHypothesesCount(void) = 0;

  constexpr static size_t max_birth_hypotheses_ = 100u;
  std::vector<ModelHypothesis> birth_hypotheses_{max_birth_hypotheses_};
};
}  // namespace mot

#endif  //  MOT_INCLUDE_BIRTH_MODEL_HPP_
