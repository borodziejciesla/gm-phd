#ifndef INCLUDE_BIRTH_GENERATOR_FACTORY_HPP_
#define INCLUDE_BIRTH_GENERATOR_FACTORY_HPP_

#include <memory>
#include <stdexcept>

#include "birth_generator_random.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class BirthGeneratorFactory {
 public:
  static std::unique_ptr<BirthGeneratorInterface<state_size, measurement_size>> Create(
      const BirthType birth_type, const BirthGeneratorCalibration& calibrations) {
    switch (birth_type) {
      case BirthType::Random:
        return std::make_unique<BirthGeneratorRandom<state_size, measurement_size>>(calibrations);
      default:
        throw std::invalid_argument("BirthGeneratorFactory::Create(): Invalid BirthType!");
    }
  }
};
}  //  namespace mot

#endif  //  INCLUDE_BIRTH_GENERATOR_FACTORY_HPP_
