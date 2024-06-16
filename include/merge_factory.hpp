#ifndef INCLUDE_BIRTH_MERGE_FACTORY_HPP_
#define INCLUDE_BIRTH_MERGE_FACTORY_HPP_

#include <memory>
#include <stdexcept>

#include "merge_classic.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class MergeFactory {
 public:
  static std::unique_ptr<MergeInterface<state_size, measurement_size>> Create(
      const MergeType merge_type) {
    switch (merge_type) {
      case MergeType::Classic:
        return std::make_unique<MergeClassic<state_size, measurement_size>>();
      default:
        throw std::invalid_argument("MergeFactory::Create(): Invalid MergeType!");
    }
  }
};
}  //  namespace mot

#endif  //  INCLUDE_BIRTH_MERGE_FACTORY_HPP_
