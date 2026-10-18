#include "gm_phd.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "birth_model_dummy.hpp"
#include "cv_motion_model.hpp"
#include "gtest/gtest.h"
#include "ospa.hpp"

namespace {
using TestMeasurement = mot::ValueWithCovariance<mot::CvMotionModel::measurement_size>;
using TestObject = mot::ValueWithCovariance<mot::CvMotionModel::state_size>;

class MeasurementBirthModel : public mot::BirthModelBase<mot::CvMotionModel> {
 protected:
  ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement) override {
    ModelHypothesis hypothesis;
    hypothesis.weight = 0.25;
    hypothesis.state(0) = measurement.value(0);
    hypothesis.state(2) = measurement.value(1);
    hypothesis.covariance = StateMatrix::Identity();
    return hypothesis;
  }

  ModelHypothesis PredictBirthHypothesis() override { return {}; }

  size_t PredictBirthHypothesesCount() override { return 0u; }
};

TestMeasurement MakeMeasurement(const double x, const double y) {
  TestMeasurement measurement;
  measurement.value << x, y;
  return measurement;
}

bool ObjectsMatchMeasurements(const std::vector<TestObject>& objects,
                              const std::vector<TestMeasurement>& measurements,
                              const double tolerance) {
  if (objects.size() != measurements.size()) {
    return false;
  }

  const auto distance = [](const TestObject& object, const TestMeasurement& measurement) {
    return std::hypot(object.value(0) - measurement.value(0),
                      object.value(2) - measurement.value(1));
  };

  std::vector<bool> matched(objects.size(), false);
  for (const auto& measurement : measurements) {
    auto closest_object = objects.size();
    auto closest_distance = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < objects.size(); ++i) {
      if (!matched[i]) {
        const auto current_distance = distance(objects[i], measurement);
        if (current_distance < closest_distance) {
          closest_object = i;
          closest_distance = current_distance;
        }
      }
    }

    if (closest_distance > tolerance) {
      return false;
    }
    matched[closest_object] = true;
  }

  return true;
}
}  // namespace

class GmPhdTests : public ::testing::Test {};

TEST_F(GmPhdTests, EvaluatesNormalizedGaussianDensity) {
  const auto density = mot::NormPdf<2u>(mot::Vector<2u>::Zero(), mot::Vector<2u>::Zero(),
                                        mot::Matrix<2u>::Identity());

  EXPECT_NEAR(density, 1.0 / (2.0 * std::numbers::pi), 1e-12);
}

TEST_F(GmPhdTests, OspaReturnsZeroForEqualPointSetsRegardlessOfOrder) {
  const std::vector<mot::OspaPosition> estimates = {{0.0, 0.0}, {5.0, 2.0}};
  const std::vector<mot::OspaPosition> truth = {{5.0, 2.0}, {0.0, 0.0}};

  EXPECT_DOUBLE_EQ(mot::OspaDistance(estimates, truth, 10.0, 2.0), 0.0);
}

TEST_F(GmPhdTests, OspaUsesMinimumCostAssignment) {
  const std::vector<mot::OspaPosition> estimates = {{0.0, 0.0}, {4.0, 0.0}};
  const std::vector<mot::OspaPosition> truth = {{3.0, 0.0}, {10.0, 0.0}};

  EXPECT_DOUBLE_EQ(mot::OspaDistance(estimates, truth, 20.0, 1.0), 4.5);
}

TEST_F(GmPhdTests, OspaIncludesLocalizationAndCardinalityError) {
  const std::vector<mot::OspaPosition> estimates = {{0.0, 0.0}, {10.0, 0.0}};
  const std::vector<mot::OspaPosition> truth = {{3.0, 0.0}};

  EXPECT_NEAR(mot::OspaDistance(estimates, truth, 10.0, 2.0), std::sqrt(54.5), 1e-12);
}

TEST_F(GmPhdTests, OspaCapsLargeLocalizationErrorsAndHandlesEmptySets) {
  const std::vector<mot::OspaPosition> estimates = {{100.0, 100.0}};
  const std::vector<mot::OspaPosition> truth = {{0.0, 0.0}};
  const std::vector<mot::OspaPosition> empty;

  EXPECT_DOUBLE_EQ(mot::OspaDistance(estimates, truth, 10.0, 2.0), 10.0);
  EXPECT_DOUBLE_EQ(mot::OspaDistance(estimates, empty, 10.0, 2.0), 10.0);
  EXPECT_DOUBLE_EQ(mot::OspaDistance(empty, empty, 10.0, 2.0), 0.0);
}

TEST_F(GmPhdTests, ConstructorTest) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, mot::BirthModelDummy>();

  EXPECT_EQ(gm_phd_filter.GetCalibration<double>("pd").value(), 0.8);
}

TEST_F(GmPhdTests, NoInputTest) {
  class BirthModelZero : public mot::BirthModelBase<mot::CvMotionModel> {
   protected:
    ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement) {
      ModelHypothesis birth_hypothesis;
      birth_hypothesis.weight = 0.25;
      birth_hypothesis.state(0) = measurement.value(0);
      birth_hypothesis.state(1) = measurement.value(1);
      birth_hypothesis.covariance = StateMatrix::Identity();

      return birth_hypothesis;
    }

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

TEST_F(GmPhdTests, CreatesSeparateObjectsForDistinctMeasurements) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();

  mot::ValueWithCovariance<2u> first_measurement;
  first_measurement.value << 0.0, 0.0;
  mot::ValueWithCovariance<2u> second_measurement;
  second_measurement.value << 10.0, 0.0;

  gm_phd_filter.Run(0.1, {first_measurement, second_measurement});

  const auto& objects = gm_phd_filter.GetObjects();
  ASSERT_EQ(objects.size(), 2u);
  EXPECT_NEAR(objects[0].value(0), 0.0, 1e-9);
  EXPECT_NEAR(objects[1].value(0), 10.0, 1e-9);
}

TEST_F(GmPhdTests, ExtractionThresholdCanBeConfigured) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();
  const auto calibration = gm_phd_filter.SetCalibrations("extraction_threshold", 0.0);
  ASSERT_TRUE(calibration.has_value());

  gm_phd_filter.Run(0.1, {MakeMeasurement(1.0, 1.0)});

  EXPECT_EQ(gm_phd_filter.GetObjects().size(), 1u);
}

TEST_F(GmPhdTests, DoesNotReuseMeasurementUpdatesWithinOneScan) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();
  const auto first_measurement = MakeMeasurement(0.0, 0.0);
  const auto second_measurement = MakeMeasurement(1.0, 1.0);

  gm_phd_filter.Run(0.1, {first_measurement, second_measurement});

  EXPECT_TRUE(std::isfinite(gm_phd_filter.GetWeightsSum()));
}

TEST_F(GmPhdTests, TracksMultipleRepeatedDetections) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();
  const auto first_measurement = MakeMeasurement(-10.0, -5.0);
  const auto second_measurement = MakeMeasurement(10.0, 5.0);

  for (int frame = 0; frame < 10; ++frame) {
    gm_phd_filter.Run(0.1 * (frame + 1), {first_measurement, second_measurement});

    const auto& objects = gm_phd_filter.GetObjects();
    ASSERT_EQ(objects.size(), 2u) << "frame " << frame;
    EXPECT_TRUE(ObjectsMatchMeasurements(objects, {first_measurement, second_measurement}, 1.0))
        << "frame " << frame;
  }
}

TEST_F(GmPhdTests, TracksFourRepeatedDetections) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();
  const std::vector<TestMeasurement> measurements = {
      MakeMeasurement(-15.0, -10.0), MakeMeasurement(-5.0, 10.0),
      MakeMeasurement(5.0, -10.0),  MakeMeasurement(15.0, 10.0)};

  for (int frame = 0; frame < 10; ++frame) {
    gm_phd_filter.Run(0.1 * (frame + 1), measurements);

    const auto& objects = gm_phd_filter.GetObjects();
    ASSERT_EQ(objects.size(), measurements.size()) << "frame " << frame;
    EXPECT_TRUE(ObjectsMatchMeasurements(objects, measurements, 1.0)) << "frame " << frame;
  }
}

TEST_F(GmPhdTests, TracksMovingTargetsInSimulation) {
  auto gm_phd_filter = mot::GmPhd<mot::CvMotionModel, MeasurementBirthModel>();

  for (int frame = 0; frame < 12; ++frame) {
    const auto first_measurement =
        MakeMeasurement(-12.0 + 0.5 * frame, -6.0 + 0.25 * frame);
    const auto second_measurement =
        MakeMeasurement(12.0 - 0.5 * frame, 6.0 - 0.25 * frame);

    gm_phd_filter.Run(0.1 * (frame + 1), {first_measurement, second_measurement});

    const auto& objects = gm_phd_filter.GetObjects();
    ASSERT_EQ(objects.size(), 2u) << "frame " << frame;
    EXPECT_TRUE(ObjectsMatchMeasurements(objects, {first_measurement, second_measurement}, 2.0))
        << "frame " << frame;
  }
}
