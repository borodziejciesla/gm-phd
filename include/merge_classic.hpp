#ifndef INCLUDE_OBJECT_MERGE_CLASSIC_HPP_
#define INCLUDE_OBJECT_MERGE_CLASSIC_HPP_

#include "merge_interface.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class MergeClassic : public MergeInterface<state_size, measurement_size> {
 public:
  using StateHypothesis = Hypothesis<state_size, measurement_size>;
  using Object = ValueWithCovariance<state_size>;

  using HypothesisVector = std::vector<StateHypothesis>;
  using ObjectsVector = std::vector<Object>;

  void MergeHypothesis(HypothesisVector& hypothesis) override {
    //
  }
};
}  // namespace mot

#endif  //  INCLUDE_OBJECT_MERGE_CLASSIC_HPP_
