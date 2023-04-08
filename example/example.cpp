#include <random>
#include <vector>

#include <iostream>

#include "matplotlibcpp.hpp"

#include "ghm_phd_cv_pose.hpp"

namespace plt = matplotlibcpp;

constexpr size_t measurements_number = 100u;
constexpr size_t objects_number = 2u;

constexpr double dt = 0.1;

std::random_device r;
std::default_random_engine e1(r());

std::uniform_real_distribution<double> pose_dist(-10.0, 10.0);
std::uniform_real_distribution<double> velocity_dist(-1.0, 1.0);

std::normal_distribution<> measurement_noisse(0.0, 0.1);

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
      measurement.covariance = 0.1 * mot::GmPhdCvPose::MeasurementSizeMatrix::Identity();
    }
  }
  return measurements;
}

/* Example */
int main () {
  //Create Filter
  mot::GmPhdCalibrations<4u, 2u> calibrations;
  
  calibrations.process_noise_diagonal = {1.0, 1.0, 100.0, 100.0};

  calibrations.observation_matrix = Eigen::Matrix<double, 2u, 4u>::Zero();
  calibrations.observation_matrix(0u, 0u) = 1.0;
  calibrations.observation_matrix(1u, 1u) = 1.0;

  calibrations.measurement_covariance = 0.2 * Eigen::Matrix<double, 2u, 2u>::Identity();

  mot::GmPhdCvPose gm_phd_filter = mot::GmPhdCvPose(calibrations);

  // Generate trajectory and measurements
  const auto trajectory = GenerateTrajectory();
  const auto measurements = GenerateMeasurements(trajectory);

  // Run Filter
  std::vector<std::vector<mot::GmPhdCvPose::Object>> objects_all;
  for (auto index = 0u; index < measurements_number; index++) {
    gm_phd_filter.Run(static_cast<double>(index) * dt, measurements.at(index));
    const auto objects = gm_phd_filter.GetObjects();
    objects_all.push_back(objects);
    std::cout << "Scan index: " << index << ", Objects number: " << objects.size() << ", Weights Sum: " << gm_phd_filter.GetWeightsSum() << "\n";
  }

  //*********************************************************************//
  // Plot
  std::vector<double> objects_x;
  std::vector<double> objects_y;
  for (auto index = 0u; index < measurements_number; index++) {
    for (const auto object : objects_all.at(index)) {
      objects_x.push_back(object.value(0u));
      objects_y.push_back(object.value(1u));
    }
  }

  std::vector<double> ref_objects_x;
  std::vector<double> ref_objects_y;
  for (auto index = 0u; index < measurements_number; index++) {
    for (const auto trajectories : trajectory.at(index)) {
      ref_objects_x.push_back(trajectories.value(0u));
      ref_objects_y.push_back(trajectories.value(1u));
    }
  }

  std::vector<double> meas_objects_x;
  std::vector<double> meas_objects_y;
  for (auto index = 0u; index < measurements_number; index++) {
    for (const auto measurement : measurements.at(index)) {
      meas_objects_x.push_back(measurement.value(0u));
      meas_objects_y.push_back(measurement.value(1u));
    }
  }

  plt::figure_size(1200, 780);
  plt::xlabel("X [m]");
  plt::ylabel("Y [m]");
  plt::plot(objects_x, objects_y, "b+");
  plt::plot(ref_objects_x, ref_objects_y, "r.");
  plt::plot(meas_objects_x, meas_objects_y, "g.");
  plt::show();

  return 0;
}
