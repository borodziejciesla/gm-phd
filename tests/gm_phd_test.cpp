#include "gm_phd.hpp"

#include <cmath>
#include <random>
#include <vector>

#include "birth_model_dummy.hpp"
#include "cv_motion_model.hpp"
#include "gtest/gtest.h"

std::mt19937 gen(0u);
std::uniform_real_distribution<> dist(-50.0, 50.0);

class GmPhdTests : public ::testing::Test {
 protected:
  void SetUp() override {
    measurements_.clear();
    timestamps_.clear();
  }

  using Measurement = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>::Measurement;
  using MeasurementsScan = std::vector<Measurement>;
  std::vector<MeasurementsScan> measurements_;
  std::vector<double> timestamps_;

  void OneMovingObjectFilterStationaryGenerateInputs() {
    // Generate measurements for a single object moving in a straight line with
    // constant velocity in the sensor frame. The object starts at (1, 1) and moves with a velocity
    // of (1, 1) m/s. The measurements are generated at a frequency of 10 Hz, and the
    // measurement noise is modeled as a Gaussian with a covariance of [[1, 0], [0, 1]] m^2. The
    // generated measurements are stored in the measurements_ member variable, which can be used for
    // testing the GmPhd filter.
    constexpr auto velocity_x = 1.0;                    // [m/s]
    constexpr auto velocity_y = 1.0;                    // [m/s]
    constexpr auto measurement_noise_covariance = 1.0;  // [m^2]
    constexpr auto dt = 0.1;                            // [s]
    constexpr auto x_start = 1.0;                       // [m]
    constexpr auto y_start = 1.0;                       // [m]
    constexpr auto num_measurements = 100u;             // [-] Number of measurements to generate

    for (auto i = 0u; i < num_measurements; ++i) {
      Measurement measurement;
      measurement.value << x_start + velocity_x * dt * i, y_start + velocity_y * dt * i;
      measurement.covariance << measurement_noise_covariance, 0.0, 0.0,
          measurement_noise_covariance;
      measurements_.push_back({measurement});
      timestamps_.push_back(dt * i);
    }
  }

  void OneMovingObjectFilterStationaryGenerateInputsWithNoise() {
    OneMovingObjectFilterStationaryGenerateInputs();

    for (auto& scan : measurements_) {
      for (auto index = 0u; index < 10u; index++) {
        Measurement measurement;
        measurement.value << dist(gen), dist(gen);
        measurement.covariance << 1.0, 0.0, 0.0, 1.0;
        scan.emplace_back(std::move(measurement));
      }
    }
  }
};

TEST_F(GmPhdTests, ConstructorTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>();

  EXPECT_EQ(gm_phd_filter.GetCalibration<double>("pd").value(), 0.8);
}

TEST_F(GmPhdTests, NoInputTest) {
  class BirthModelZero : public mot::BirthModelBase<4u, 2u> {
   protected:
    ModelHypothesis PredictBirthHypothesis(void) {
      ModelHypothesis birth_hypothesis;
      birth_hypothesis.weight = 0.0;
      birth_hypothesis.state = StateVector::Zero();
      birth_hypothesis.covariance = StateMatrix::Identity();

      return birth_hypothesis;
    }

    size_t PredictBirthHypothesesCount(void) { return 1u; }
  };

  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, BirthModelZero>();
  gm_phd_filter.Run(0.1, {});

  EXPECT_TRUE(gm_phd_filter.GetObjects().empty());
  EXPECT_DOUBLE_EQ(gm_phd_filter.GetWeightsSum(), 0.0);
}

TEST_F(GmPhdTests, AddSameInputTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>();

  mot::ValueWithCovariance<2u> measurement;
  measurement.value << 1.0, 1.0;
  measurement.covariance << 1.0, 0.0, 0.0, 1.0;

  for (int i = 0; i < 50; ++i) {
    gm_phd_filter.Run(0.1 * (i + 1), {measurement});
  }

  EXPECT_EQ(gm_phd_filter.GetObjects().size(), 1u);
  EXPECT_DOUBLE_EQ(std::round(gm_phd_filter.GetWeightsSum()), 1.0);
}

TEST_F(GmPhdTests, OneMovingObjectFilterStationaryTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>();

  OneMovingObjectFilterStationaryGenerateInputs();

  for (auto i = 0u; i < measurements_.size(); ++i) {
    gm_phd_filter.Run(timestamps_[i], measurements_[i]);

    if (i > 4u) {
      EXPECT_EQ(gm_phd_filter.GetObjects().size(), 1u)
          << std::string("Failed at index: ") << std::to_string(i);
      EXPECT_DOUBLE_EQ(std::round(gm_phd_filter.GetWeightsSum()), 1.0);

      const auto& object = gm_phd_filter.GetObjects().front();
      EXPECT_NEAR(object.value(0), measurements_[i].front().value(0), 0.2);
      EXPECT_NEAR(object.value(1), 1.0, 0.2);
      EXPECT_NEAR(object.value(2), measurements_[i].front().value(1), 0.2);
      EXPECT_NEAR(object.value(3), 1.0, 0.2);
    }
  }
}

TEST_F(GmPhdTests, OneMovingObjectFilterStationaryTestWithNoise) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>();

  OneMovingObjectFilterStationaryGenerateInputsWithNoise();

  auto weights_sum = 0.0;

  for (auto i = 0u; i < measurements_.size(); ++i) {
    gm_phd_filter.Run(timestamps_[i], measurements_[i]);
    weights_sum += gm_phd_filter.GetWeightsSum();

    if (gm_phd_filter.GetObjects().size() > 1u) {
      const auto& object = gm_phd_filter.GetObjects().front();
      EXPECT_NEAR(object.value(0), measurements_[i].front().value(0), 0.1);
      EXPECT_NEAR(object.value(2), measurements_[i].front().value(1), 0.1);
    }
  }

  EXPECT_DOUBLE_EQ(std::round(weights_sum / measurements_.size()), 1.0);
}