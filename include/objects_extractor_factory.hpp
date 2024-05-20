#ifndef INCLUDE_OBJECT_EXTRACTOR_FACTORY_HPP_
#define INCLUDE_OBJECT_EXTRACTOR_FACTORY_HPP_

#include <memory>

#include "object_extractor_classic.hpp"
#include "object_extractor_interface.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class ObjectsExtractorFactory {
 public:
  static std::unique_ptr<ObjectExtractorInterface<state_size, measurement_size>>
  CreateExtractor() {  // TODO: Add calibration
    return std::make_unique<ObjectExtractorClassic<state_size, measurement_size>>();
  }
};
}  //  namespace mot

#endif  //  INCLUDE_OBJECT_EXTRACTOR_FACTORY_HPP_
