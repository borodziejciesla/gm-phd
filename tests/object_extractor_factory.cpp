#include <cmath>
#include <memory>
#include <numbers>

#include "../include/aliases.hpp"
#include "../include/objects_extractor_factory.hpp"
#include "gtest/gtest.h"

class ObjectsExtractorTests : public ::testing::Test {
 protected:
  void SetUp(void) override {}
};

TEST_F(ObjectsExtractorTests, GenerateSuccessRandomTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  using Factory = mot::ObjectsExtractorFactory<state_size, measurement_size>;
  EXPECT_NO_THROW(Factory::Create(mot::ExtractorType::Classic));
}

TEST_F(ObjectsExtractorTests, GenerateFailRandomTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  using Factory = mot::ObjectsExtractorFactory<state_size, measurement_size>;
  EXPECT_THROW(Factory::Create(static_cast<mot::ExtractorType>(1u)), std::invalid_argument);
}
