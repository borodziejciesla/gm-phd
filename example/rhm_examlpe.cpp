#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <string>

#include "matplotlibcpp.hpp"

#include "rhm_gm_phd.hpp"
#include "trajectory_generation.hpp"

namespace plt = matplotlibcpp;

/*************************** Main ***************************/
int main() {
  const auto gt = GetGroundTruth();

  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(5.0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rand(0.0, 1.0);

  std::normal_distribution<double> normal_distribution(0.0, 1.0);

  /************************** Define tracker object **************************/
  mot::GmPhdCalibrations<4u, 2u> calibrations;
  calibrations.process_noise_diagonal = {1.0, 1.0, 100.0, 100.0};
  calibrations.observation_matrix = Eigen::Matrix<double, 2u, 4u>::Zero();
  calibrations.observation_matrix(0u, 0u) = 1.0;
  calibrations.observation_matrix(1u, 1u) = 1.0;
  calibrations.measurement_covariance = 0.2 * Eigen::Matrix<double, 2u, 2u>::Identity();

  mot::RhmGmPhd gm_phd_filter = mot::RhmGmPhd();

  /************************** Run **************************/
  using Measurements = std::vector<mot::RhmGmPhd::Measurement>;

  std::vector<std::vector<mot::RhmGmPhd::Object>> output_objects;
  std::vector<std::vector<mot::RhmGmPhd::Measurement>> detections;

  for (auto index = 0u; index < gt.time_steps; index++) {
    // Select detctions number in step
    auto detections_number = distribution(generator);
    while (detections_number == 0)
      detections_number = distribution(generator);

    std::cout << "Time step: " << std::to_string(index) << ", " << std::to_string(detections_number) << " Measurements\n";

    // Generate noisy measurement
    std::vector<mot::RhmGmPhd::Measurement> measurements(detections_number);
    for (auto & measurement : measurements) {
      std::array<double, 2u> h = {-1.0 + 2 * rand(gen), -1.0 + 2 * rand(gen)};
      while (std::hypot(h.at(0), h.at(1)) > 1.0)
        h = {-1.0 + 2 * rand(gen), -1.0 + 2 * rand(gen)};

      measurement.value(0u) = gt.center.at(index).at(0u)
        + h.at(0) * gt.size.at(index).first * std::cos(gt.orientation.at(index))
        - h.at(1) * gt.size.at(index).second * std::sin(gt.orientation.at(index));// + 10.0 * normal_distribution(generator);
      measurement.value(1u) = gt.center.at(index).at(1u)
        + h.at(0) * gt.size.at(index).first * std::sin(gt.orientation.at(index))
        + h.at(1) * gt.size.at(index).second * std::cos(gt.orientation.at(index));// + 10.0 * normal_distribution(generator);

      measurement.covariance = Eigen::Matrix<double, mot::measurement_size, mot::measurement_size>::Zero();
      measurement.covariance(0u, 0u) = 200.0;
      measurement.covariance(1u, 1u) = 8.0;
    }

    detections.push_back(measurements);

    // Run algo
    gm_phd_filter.Run(static_cast<double>(index) * 10.0, measurements);
    const auto objects_output = gm_phd_filter.GetObjects();

    if (!objects_output.empty()) {
        output_objects.push_back(objects_output);
        std::cout << "Objects number: " << objects_output.size() << "\n";
    }
  }

  /************************** Plot outputs **************************/
  plt::figure_size(1200, 780);

  plt::xlabel("X [m]");
  plt::ylabel("Y [m]");

  // Trajectory
  std::vector<double> x_traj;
  std::vector<double> y_traj;
  for (const auto & point : gt.center) {
    x_traj.push_back(point.at(0u));
    y_traj.push_back(point.at(1u));
  }

  std::map<std::string, std::string> keywords_traj;
  keywords_traj.insert(std::pair<std::string, std::string>("label", "Trajectory") );

  plt::plot(x_traj, y_traj, keywords_traj);

  // Detections
  std::vector<double> x_detections;
  std::vector<double> y_detections;
  for (auto index = 0u; index < detections.size(); index = index + 3u) {
    for (const auto & detection : detections.at(index)) {
      x_detections.push_back(detection.value(0u));
      y_detections.push_back(detection.value(1u));
    }
  }
  plt::scatter(x_detections, y_detections);

  // Objects Center
//   std::vector<double> x_objects;
//   std::vector<double> y_objects;
//   for (auto index = 0u; index < output_objects.size(); index = index + 3u) {
//     x_objects.push_back(output_objects.at(index).kinematic_state.value(0u));
//     y_objects.push_back(output_objects.at(index).kinematic_state.value(1u));

//     const auto [x_ellips, y_ellipse] = CreateEllipse(output_objects.at(index).extent_state.extent_state, std::make_pair(output_objects.at(index).kinematic_state.state(0u), output_objects.at(index).kinematic_state.state(1u)));
//     plt::plot(x_ellips, y_ellipse, "r");
//   }
//   plt::plot(x_objects, y_objects, "r*");

//   // Reference
//   std::vector<double> x_ref;
//   std::vector<double> y_ref;
//   for (auto index = 0u; index < gt.time_steps; index = index + 3u) {
//     x_ref.push_back(gt.center.at(index).at(0u));
//     y_ref.push_back(gt.center.at(index).at(1u));

//     mot::ExtentState ellipse = {gt.orientation.at(index), gt.size.at(index).first, gt.size.at(index).second};
//     const auto [x_ellips, y_ellipse] = CreateEllipse(ellipse, std::make_pair(gt.center.at(index).at(0u), gt.center.at(index).at(1u)));
//     plt::plot(x_ellips, y_ellipse, "k");
//   }
//   plt::plot(x_ref, y_ref, "k*");
//   plt::show();

//   // Velocity
//   std::vector<double> vx_ref;
//   std::vector<double> vy_ref;
//   std::vector<double> vx_obj;
//   std::vector<double> vy_obj;
//   std::vector<double> idx;
//   for (auto index = 0u; index < gt.time_steps; index = index + 1u) {
//     vx_ref.push_back(gt.velocity.at(index).first);
//     vy_ref.push_back(gt.velocity.at(index).second);

//     vx_obj.push_back(output_objects.at(index).kinematic_state.state(2u));
//     vy_obj.push_back(output_objects.at(index).kinematic_state.state(3u));

//     idx.push_back(index);
//   }
//   plt::plot(idx, vx_ref, "b");
//   plt::plot(idx, vy_ref, "b:");
//   plt::plot(idx, vx_obj, "r");
//   plt::plot(idx, vy_obj, "r:");
//   plt::show();

  return EXIT_SUCCESS;
}