#ifndef INCLUDE_BIRTH_GENERATOR_INTERFACE_HPP_
#define INCLUDE_BIRTH_GENERATOR_INTERFACE_HPP_

#include <vector>

#include "gm_phd_calibrations.hpp"
#include "hypothesis.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class BirthGeneratorInterface {
 public:
  explicit BirthGeneratorInterface(const BirthGeneratorCalibration& calibrations) {}

  virtual const std::vector<Hypothesis<state_size, measurement_size>>& Run(void) = 0;
};
}  //  namespace mot

#endif  //  INCLUDE_BIRTH_GENERATOR_INTERFACE_HPP_
