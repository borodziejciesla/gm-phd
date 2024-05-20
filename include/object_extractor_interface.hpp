#ifndef INCLUDE_OBJECT_EXTRACTOR_HPP_
#define INCLUDE_OBJECT_EXTRACTOR_HPP_

#include <vector>

#include "aliases.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class ObjectExtractorInterface {
 public:
  using StateHypothesis = Hypothesis<state_size, measurement_size>;
  using Object = ValueWithCovariance<state_size>;

  using HypothesisVector = std::vector<StateHypothesis>;
  using ObjectsVector = std::vector<Object>;

  virtual const ObjectsVector& ExtractObjects(const HypothesisVector& hypothesis) = 0;
};
}  // namespace mot

#endif  //  INCLUDE_OBJECT_EXTRACTOR_HPP_
