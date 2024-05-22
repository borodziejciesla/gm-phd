#include "../include/gm_phd.hpp"

#include <memory>
#include <numbers>

#include "gtest/gtest.h"

class GmPhdTests : public ::testing::Test {
 protected:
  void SetUp(void) override {}
};

TEST_F(GmPhdTests, ConstructorTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::GmPhdCalibrations<state_size, measurement_size> calibrations;

  auto gm_phd = std::make_unique<mot::GmPhd<state_size, measurement_size>>(calibrations);
}
