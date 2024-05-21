#include "../include/hypothesis.hpp"

#include <memory>
#include <numbers>

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

TEST_F(HypothesisTests, PredictTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  auto predict_hypothesis = [](mot::Hypothesis<state_size, measurement_size>& hypothesis,
                               float time_delta) {
    hypothesis.State() = time_delta * hypothesis.State();
    hypothesis.Covariance() = time_delta * hypothesis.Covariance();
    hypothesis.Weight() = 1.0f;
  };

  constexpr auto weight = 0.66f;
  const auto state = mot::Vector<state_size>::Ones();
  const auto covariance = mot::SquareMatrix<state_size>::Identity();

  auto hypothesis =
      std::make_unique<mot::Hypothesis<state_size, measurement_size>>(weight, state, covariance);
  hypothesis->PredictState(predict_hypothesis, 2.0f);

  EXPECT_FLOAT_EQ(hypothesis->Weight(), 1.0f);
  EXPECT_EQ(hypothesis->State(), state * 2.0f);
  EXPECT_EQ(hypothesis->Covariance(), covariance * 2.0f);
}

TEST_F(HypothesisTests, MoveSensorTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  constexpr auto weight = 0.66f;
  const auto state = mot::Vector<state_size>::Ones();
  const auto covariance = mot::SquareMatrix<state_size>::Ones();

  const mot::SensorPoseVector sensor_pose_delta = {2.0f, 3.0f, std::numbers::pi_v<float>};

  auto hypothesis =
      std::make_unique<mot::Hypothesis<state_size, measurement_size>>(weight, state, covariance);
  hypothesis->MoveSensor(sensor_pose_delta);

  const mot::Vector<state_size> result = {1.0f, 2.0f, 1.0f, 1.0f};

  EXPECT_FLOAT_EQ(hypothesis->Weight(), weight);
  EXPECT_FLOAT_EQ(hypothesis->State()[0], result[0]);
  EXPECT_FLOAT_EQ(hypothesis->State()[1], result[1]);
  EXPECT_EQ(hypothesis->Covariance(), covariance);
}
