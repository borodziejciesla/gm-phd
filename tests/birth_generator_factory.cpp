#include "../include/birth_generator_factory.hpp"

#include <cmath>
#include <memory>
#include <numbers>

#include "../include/aliases.hpp"
#include "gtest/gtest.h"

class BirthGeneratorTests : public ::testing::Test {
 protected:
  void SetUp(void) override {}
};

TEST_F(BirthGeneratorTests, GenerateSuccessRandomTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::BirthGeneratorCalibration calibrations;
  mot::BirthType birth_type = mot::BirthType::Random;

  using Factory = mot::BirthGeneratorFactory<state_size, measurement_size>;
  EXPECT_NO_THROW(Factory::Create(birth_type, calibrations));
}

TEST_F(BirthGeneratorTests, GenerateFailRandomTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::BirthGeneratorCalibration calibrations;
  mot::BirthType birth_type = static_cast<mot::BirthType>(1u);

  using Factory = mot::BirthGeneratorFactory<state_size, measurement_size>;
  EXPECT_THROW(Factory::Create(birth_type, calibrations), std::invalid_argument);
}
