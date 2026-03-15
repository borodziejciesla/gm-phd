#include "gm_phd.hpp"

#include "cv_motion_model.hpp"
#include "gtest/gtest.h"

class GmPhdTests : public ::testing::Test {};

TEST_F(GmPhdTests, ConstructorTest) {
  constexpr auto state_size = 4u;
  constexpr auto observation_size = 2u;

  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel>();
}
