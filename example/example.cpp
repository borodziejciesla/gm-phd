#include <random>
#include <vector>

#include "ghm_phd_cv_pose.hpp"

constexpr size_t measurements_number = 100u;
constexpr size_t objects_number = 2u;

constexpr double dt = 0.1;

std::random_device r;
std::default_random_engine e1(r());

std::uniform_real_distribution<double> pose_dist(-150, 150);
std::uniform_real_distribution<double> velocity_dist(-10, 10);

std::normal_distribution<> measurement_noisse(0.0, 1.0);

using Measurements = std::vector<mot::GmPhdCvPose::Measurement>;

std::vector<Measurements> GenerateTrajectory(void) {
  std::vector<Measurements> trajectory(measurements_number);

  // Initialize
  std::vector<mot::GmPhdCvPose::Object> objects(objects_number);
  for (auto & object : objects) {
    object.value(0u) = pose_dist(e1);
    object.value(1u) = pose_dist(e1);
    object.value(2u) = velocity_dist(e1);
    object.value(3u) = velocity_dist(e1);
  }

  // Generate
  for (auto index = 0u; index < measurements_number; index++) {
    const auto time = static_cast<double>(index) * dt;
    for (auto object_index = 0u; object_index < objects_number; object_index++) {
      mot::GmPhdCvPose::Measurement measurement;
      measurement.value(0u) = objects.at(object_index).value(0u) + time * objects.at(object_index).value(2u);
      measurement.value(1u) = objects.at(object_index).value(1u) + time * objects.at(object_index).value(3u);

      trajectory.at(index).push_back(measurement);
    }
  }

  return trajectory;
}

std::vector<Measurements> GenerateMeasurements(const std::vector<Measurements> & trajectory) {
  std::vector<Measurements> measurements = trajectory;
  for (auto & measurement_scan : measurements) {
    for (auto & measurement : measurement_scan) {
      measurement.value(0u) += measurement_noisse(e1);
      measurement.value(1u) += measurement_noisse(e1);
      measurement.covariance = mot::GmPhdCvPose::MeasurementSizeMatrix::Identity();
    }
  }
  return measurements;
}

/* Example */
int main () {
  //Create Filter
  mot::GmPhdCalibrations<4u, 2u> calibrations;
  mot::GmPhdCvPose gm_phd_filter = mot::GmPhdCvPose(calibrations);

  // Generate trajectory and measurements
  const auto trajectory = GenerateTrajectory();
  const auto measurements = GenerateMeasurements(trajectory);

  // Run Filter
  for (auto index = 0u; index < measurements_number; index++) {
    gm_phd_filter.Run(static_cast<double>(index) * dt, measurements.at(index));
    const auto objects = gm_phd_filter.GetObjects();
  }

  return 0;
}
