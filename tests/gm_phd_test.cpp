#include "gm_phd.hpp"

#include <cmath>

#include "cv_motion_model.hpp"
#include "gtest/gtest.h"

class GmPhdTests : public ::testing::Test {};

TEST_F(GmPhdTests, ConstructorTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel>();

  EXPECT_EQ(gm_phd_filter.GetCalibration<double>("pd").value(), 0.8);
}

TEST_F(GmPhdTests, NoInputTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel>();
  gm_phd_filter.Run(0.1, {});

  EXPECT_TRUE(gm_phd_filter.GetObjects().empty());
  EXPECT_DOUBLE_EQ(gm_phd_filter.GetWeightsSum(), 0.0);
}

TEST_F(GmPhdTests, AddSameInputTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel>();

  mot::ValueWithCovariance<2u> measurement;
  measurement.value << 1.0, 1.0;
  measurement.covariance << 1.0, 0.0, 0.0, 1.0;

  for (int i = 0; i < 50; ++i) {
    gm_phd_filter.Run(0.1 * (i + 1), {measurement});
  }

  EXPECT_EQ(gm_phd_filter.GetObjects().size(), 1u);
  EXPECT_DOUBLE_EQ(std::round(gm_phd_filter.GetWeightsSum()), 1.0);
}