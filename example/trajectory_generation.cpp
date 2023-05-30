#include "trajectory_generation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

GroundTruth GetGroundTruth(void) {
  // Generate orientation
  std::vector<double> gt_orient(103);
  std::fill(gt_orient.begin(), gt_orient.begin() + 20u, -std::numbers::pi_v<double>/4); // 1
  std::generate(gt_orient.begin() + 20u, gt_orient.begin() + 31u, [n = 0]() mutable { return (-std::numbers::pi_v<double>/4 + (static_cast<double>(n++) * std::numbers::pi_v<double> / 40.0)); });    // 1
  std::fill(gt_orient.begin() + 31u, gt_orient.begin() + 41u, 0.0);   // 3
  std::generate(gt_orient.begin() + 41u, gt_orient.begin() + 52u, [n = 0]() mutable { return (static_cast<double>(n++) * std::numbers::pi_v<double> / 20.0); });  // 4
  std::fill(gt_orient.begin() + 52u, gt_orient.begin() + 72u, std::numbers::pi_v<double> / 2.0);   // 5
  std::generate(gt_orient.begin() + 72u, gt_orient.begin() + 83u, [n = 0]() mutable { return (std::numbers::pi_v<double> / 2.0 + static_cast<double>(n++) * std::numbers::pi_v<double> / 20.0); });    // 6
  std::fill(gt_orient.begin() + 83u, gt_orient.end(), std::numbers::pi_v<double>);   // 7

  // Generate vel
  std::vector<std::pair<double, double>> gt_vel;
  std::transform(gt_orient.begin(), gt_orient.end(),
    std::back_inserter(gt_vel),
    [](const double orientation) {
      return std::make_pair<double, double>((500.0 / 36.0) * std::cos(orientation), (500.0 / 36.0) * std::sin(orientation));
    }
  );

  // Generate size
  std::vector<std::pair<double, double>> gt_length(103);
  std::fill(gt_length.begin(), gt_length.end(), std::make_pair<double, double>(170.0, 40.0));

  const auto time_steps = gt_orient.size();
  constexpr auto time_interval = 10.0;

  // Rotation
  std::vector<std::array<double, 4u>> gt_rotation;
  std::transform(gt_orient.begin(), gt_orient.end(),
    std::back_inserter(gt_rotation),
    [](const double orientation)  {
      std::array<double, 4u> rotation = {std::cos(orientation), -std::sin(orientation), std::sin(orientation), std::cos(orientation)};
      return rotation;
    }
  );

  // Center
  std::vector<std::array<double, 2u>> gt_center(time_steps);
  gt_center.at(0) = {0.0, 0.0};

  for (auto index = 1u; index < time_steps; index++) {
    gt_center.at(index).at(0u) = gt_center.at(index - 1u).at(0) + gt_vel.at(index).first * time_interval;
    gt_center.at(index).at(1u) = gt_center.at(index - 1u).at(1) + gt_vel.at(index).second * time_interval;
  }

  // Set output
  GroundTruth gt;
  gt.time_steps = time_steps;
  gt.orientation = gt_orient;
  gt.velocity = gt_vel;
  gt.size = gt_length;
  gt.rotation = gt_rotation;
  gt.center = gt_center;

  return gt;
}
