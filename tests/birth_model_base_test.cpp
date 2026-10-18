#include "birth_model_base.hpp"

#include "cv_motion_model.hpp"
#include "gtest/gtest.h"

namespace {
class TestBirthModel : public mot::BirthModelBase<mot::CvMotionModel> {
 public:
  void SetPredictedHypothesesCount(const size_t count) { predicted_hypotheses_count_ = count; }

 protected:
  ModelHypothesis PredictBirthHypothesisFromMeasurement(const Measurement& measurement) override {
    ModelHypothesis hypothesis;
    hypothesis.weight = measurement.value(0);
    hypothesis.state(0) = measurement.value(0);
    hypothesis.state(1) = measurement.value(1);
    return hypothesis;
  }

  ModelHypothesis PredictBirthHypothesis() override {
    ModelHypothesis hypothesis;
    hypothesis.weight = 0.75;
    hypothesis.state = StateVector::Ones();
    return hypothesis;
  }

  size_t PredictBirthHypothesesCount() override { return predicted_hypotheses_count_; }

 private:
  size_t predicted_hypotheses_count_ = 2u;
};

mot::ValueWithCovariance<mot::CvMotionModel::measurement_size> MakeMeasurement(
    const double x, const double y) {
  mot::ValueWithCovariance<mot::CvMotionModel::measurement_size> measurement;
  measurement.value << x, y;
  return measurement;
}
}  // namespace

TEST(BirthModelBaseTests, CreatesHypothesesFromMeasurementsAndModelPredictions) {
  TestBirthModel birth_model;

  const auto& hypotheses = birth_model.PreparedictBirthHypotheses(
      {MakeMeasurement(1.0, 2.0), MakeMeasurement(3.0, 4.0)});

  ASSERT_EQ(hypotheses.size(), 4u);
  EXPECT_DOUBLE_EQ(hypotheses[0].weight, 1.0);
  EXPECT_DOUBLE_EQ(hypotheses[0].state(0), 1.0);
  EXPECT_DOUBLE_EQ(hypotheses[0].state(1), 2.0);
  EXPECT_DOUBLE_EQ(hypotheses[1].weight, 3.0);
  EXPECT_DOUBLE_EQ(hypotheses[1].state(0), 3.0);
  EXPECT_DOUBLE_EQ(hypotheses[1].state(1), 4.0);
  EXPECT_DOUBLE_EQ(hypotheses[2].weight, 0.75);
  EXPECT_TRUE(hypotheses[2].state.isOnes());
  EXPECT_DOUBLE_EQ(hypotheses[3].weight, 0.75);
  EXPECT_TRUE(hypotheses[3].state.isOnes());
}

TEST(BirthModelBaseTests, ReplacesPreviouslyPreparedHypotheses) {
  TestBirthModel birth_model;

  ASSERT_EQ(birth_model.PreparedictBirthHypotheses({MakeMeasurement(1.0, 2.0)}).size(), 3u);

  birth_model.SetPredictedHypothesesCount(1u);
  const auto& hypotheses = birth_model.PreparedictBirthHypotheses({});

  ASSERT_EQ(hypotheses.size(), 1u);
  EXPECT_DOUBLE_EQ(hypotheses[0].weight, 0.75);
  EXPECT_TRUE(hypotheses[0].state.isOnes());
}
