#ifndef INCLUDE_OBJECT_EXTRACTOR_CLASSIC_HPP_
#define INCLUDE_OBJECT_EXTRACTOR_CLASSIC_HPP_

#include "object_extractor_interface.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class ObjectExtractorClassic : public ObjectExtractorInterface<state_size, measurement_size> {
 public:
  using StateHypothesis = Hypothesis<state_size, measurement_size>;
  using Object = ValueWithCovariance<state_size>;

  using HypothesisVector = std::vector<StateHypothesis>;
  using ObjectsVector = std::vector<Object>;

  const ObjectsVector& ExtractObjects(const HypothesisVector& hypothesis) {
    static ObjectsVector objects;
    return objects;
  }
};
}  // namespace mot

#endif  //  INCLUDE_OBJECT_EXTRACTOR_CLASSIC_HPP_
