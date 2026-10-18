#ifndef MOT_INCLUDE_BIRTH_MODEL_HPP_
#define MOT_INCLUDE_BIRTH_MODEL_HPP_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "helpers.hpp"
#include "hypothesis.hpp"

namespace mot {
template <typename MotionModel>
class BirthModelBase {
 public:
  using Measurement = ValueWithCovariance<MotionModel::measurement_size>;

  using StateVector = Vector<MotionModel::state_size>;
  using StateMatrix = Matrix<MotionModel::state_size>;

  using ModelHypothesis = Hypothesis<MotionModel::state_size, MotionModel::measurement_size>;

  const std::vector<ModelHypothesis>& PreparedictBirthHypotheses(
      const std::vector<Measurement>& measurements) {
    birth_hypotheses_.clear();

    std::transform(measurements.begin(), measurements.end(), std::back_inserter(birth_hypotheses_),
                   [this](const Measurement& measurement) {
                     return PredictBirthHypothesisFromMeasurement(measurement);
                   });

    for (size_t i = 0u; i < PredictBirthHypothesesCount(); ++i) {
      birth_hypotheses_.emplace_back(std::move(PredictBirthHypothesis()));
    }

    return birth_hypotheses_;
  }

 protected:
  virtual ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement) = 0;
  virtual ModelHypothesis PredictBirthHypothesis(void) = 0;
  virtual size_t PredictBirthHypothesesCount(void) = 0;

  constexpr static size_t max_birth_hypotheses_ = 10u;
  std::vector<ModelHypothesis> birth_hypotheses_{max_birth_hypotheses_};
};
}  // namespace mot

#endif  //  MOT_INCLUDE_BIRTH_MODEL_HPP_
