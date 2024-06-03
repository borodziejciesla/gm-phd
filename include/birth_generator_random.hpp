#ifndef INCLUDE_BIRTH_GENERATOR_RANDOM_HPP_
#define INCLUDE_BIRTH_GENERATOR_RANDOM_HPP_

#include "birth_generator_interface.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
class BirthGeneratorRandom : public BirthGeneratorInterface<state_size, measurement_size> {
 public:
  explicit BirthGeneratorRandom(const BirthGeneratorCalibration& calibrations)
      : BirthGeneratorInterface<state_size, measurement_size>(calibrations),
        calibrations_{calibrations},
        new_hypothesis_{
            std::vector<Hypothesis<state_size, measurement_size>>(calibrations_.birth_number)} {
  }  // TODO: check issue with inheritance

  const std::vector<Hypothesis<state_size, measurement_size>>& Run(void) override final {
    for (auto index = 0u; index < calibrations_.birth_number; index++) {
      new_hypothesis_.at(index) =
          Hypothesis<state_size, measurement_size>();  // TODO: implement randomized
    }

    return new_hypothesis_;
  }

 protected:
  BirthGeneratorCalibration calibrations_;
  std::vector<Hypothesis<state_size, measurement_size>> new_hypothesis_;
};
}  //  namespace mot

#endif  //  INCLUDE_BIRTH_GENERATOR_RANDOM_HPP_
