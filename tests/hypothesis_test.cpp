#include "../include/hypothesis.hpp"

#include <memory>

#include "gtest/gtest.h"

class HypothesisTests : public ::testing::Test {
 protected:
  void SetUp(void) override {}
};

TEST_F(HypothesisTests, DefaultConstructorTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  auto hypothesis = std::make_unique<mot::Hypothesis<state_size, measurement_size>>();
  EXPECT_EQ(hypothesis->Weight(), 0.0f);
  EXPECT_EQ(hypothesis->State(), mot::Vector<state_size>::Zero());
  EXPECT_EQ(hypothesis->Covariance(), mot::SquareMatrix<state_size>::Zero());
}

TEST_F(HypothesisTests, CustomConstructorTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  constexpr auto weight = 0.66f;
  const auto state = mot::Vector<state_size>::Ones();
  const auto covariance = mot::SquareMatrix<state_size>::Ones();

  auto hypothesis =
      std::make_unique<mot::Hypothesis<state_size, measurement_size>>(weight, state, covariance);
  EXPECT_FLOAT_EQ(hypothesis->Weight(), weight);
  EXPECT_EQ(hypothesis->State(), state);
  EXPECT_EQ(hypothesis->Covariance(), covariance);
}
