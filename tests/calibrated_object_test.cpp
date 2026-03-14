#include "calibrated_object.hpp"

#include "gtest/gtest.h"

class CalibratedObjectTests : public ::testing::Test {};

TEST_F(CalibratedObjectTests, ConstructorTest) {
  std::unique_ptr<mot::CalibratedObject> co;
  EXPECT_NO_THROW(co = std::make_unique<mot::CalibratedObject>());
}

// Set only invalid - no calibration added to map
TEST_F(CalibratedObjectTests, SetCalibrationInvalidTest) {
  mot::CalibratedObject co;
  EXPECT_FALSE(co.SetCalibrations("invalid", 1.0).has_value());
  EXPECT_FALSE(co.SetCalibrations("invalid2", 1.0).has_value());
  EXPECT_FALSE(co.SetCalibrations("invalid3", 1.0).has_value());
  EXPECT_FALSE(co.SetCalibrations("invalid4", 1.0).has_value());
}

TEST_F(CalibratedObjectTests, SetMultipleCalibrationInvalidTest) {
  mot::CalibratedObject co;
  EXPECT_FALSE(
      co.SetCalibrations("invalid", 1.0, "also_invalid", 2.0, "third_invalid", 3.0).has_value());
  EXPECT_FALSE(
      co.SetCalibrations("invalid1", 1.0, "also_invalid1", 2.0, "third_invalid1", 3.0).has_value());
}

// Get only invalid - no calibration added to map
TEST_F(CalibratedObjectTests, GetCalibrationInvalidTest) {
  mot::CalibratedObject co;
  EXPECT_FALSE(co.GetCalibration<double>("invalid").has_value());
  EXPECT_FALSE(co.GetCalibration<double>("invalid2").has_value());
  EXPECT_FALSE(co.GetCalibration<double>("invalid3").has_value());
  EXPECT_FALSE(co.GetCalibration<double>("invalid4").has_value());
}
