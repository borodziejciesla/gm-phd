#include "../include/gm_phd.hpp"

#include <cmath>
#include <memory>
#include <numbers>

#include "../include/aliases.hpp"
#include "gtest/gtest.h"

class GmPhdTests : public ::testing::Test {
 protected:
  void SetUp(void) override {}
};

TEST_F(GmPhdTests, ConstructorTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::GmPhdCalibrations<state_size, measurement_size> calibrations;

  calibrations.predict_hypothesis = [](mot::Hypothesis<state_size, measurement_size>& hypothesis,
                                       const float time_delta) {
    // Do nothing
  };

  calibrations.predict_observation = [](mot::Hypothesis<state_size, measurement_size>& hypothesis) {
    // Do nothing
  };

  auto gm_phd = std::make_unique<mot::GmPhd<state_size, measurement_size>>(calibrations);
}

TEST_F(GmPhdTests, SimpleRunTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::GmPhdCalibrations<state_size, measurement_size> calibrations;
  mot::GmPhd<state_size, measurement_size>::MeasurementsVector measurements;

  calibrations.predict_hypothesis = [state_size, measurement_size](
                                        mot::Hypothesis<state_size, measurement_size>& hypothesis,
                                        const float time_delta) {
    // Do nothing
  };

  calibrations.predict_observation =
      [state_size, measurement_size](mot::Hypothesis<state_size, measurement_size>& hypothesis) {
        // Do nothing
      };

  auto gm_phd = std::make_unique<mot::GmPhd<state_size, measurement_size>>(calibrations);
  EXPECT_NO_THROW(gm_phd->Run(0.0, measurements));
}

TEST_F(GmPhdTests, CvRunTest) {
  constexpr auto state_size = 4u;
  constexpr auto measurement_size = 2u;

  mot::GmPhdCalibrations<state_size, measurement_size> calibrations;
  mot::GmPhd<state_size, measurement_size>::MeasurementsVector measurements;

  calibrations.predict_hypothesis = [](mot::Hypothesis<state_size, measurement_size>& hypothesis,
                                       const float time_delta) {
    // Set transition matrix
    static mot::SquareMatrix<4u> transition_matrix = mot::SquareMatrix<4u>::Identity();
    transition_matrix(0u, 1u) = time_delta;
    transition_matrix(2u, 3u) = time_delta;

    // Set new state
    hypothesis.State() = transition_matrix * hypothesis.State();

    // Set process noise covariance matrix
    static mot::SquareMatrix<2u> q = mot::SquareMatrix<2u>::Identity();
    q(0u, 0u) = 1.0f;
    q(1u, 1u) = 1.0f;
    static mot::Matrix<4u, 2u> g = mot::Matrix<4u, 2u>::Zero();
    g(0u, 0u) = 0.5f * static_cast<float>(std::pow(time_delta, 2u));
    g(1u, 0u) = time_delta;
    g(2u, 1u) = 0.5f * static_cast<float>(std::pow(time_delta, 2u));
    g(3u, 1u) = time_delta;
    static mot::SquareMatrix<4u> pnc = mot::SquareMatrix<4u>::Identity();
    pnc = g * q * g.transpose();

    // Set state covariance
    hypothesis.Covariance() =
        transition_matrix * hypothesis.Covariance() * transition_matrix.transpose() + pnc;
  };

  calibrations.predict_observation = [](mot::Hypothesis<state_size, measurement_size>& hypothesis) {
    // Set observation matrix
    static mot::Matrix<2u, 4u> observation_matrix = mot::Matrix<2u, 4u>::Zero();
    observation_matrix(0u, 0u) = 1.0f;
    observation_matrix(1u, 2u) = 1.0f;

    // Set observation
    hypothesis.PredictedMeasurement() = observation_matrix * hypothesis.State();
  };

  auto gm_phd = std::make_unique<mot::GmPhd<state_size, measurement_size>>(calibrations);
  EXPECT_NO_THROW(gm_phd->Run(0.0, measurements));
}
