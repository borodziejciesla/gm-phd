#ifndef GM_PHD_EXAMPLE_TRAJECTORY_GENERATION_HPP_
#define GM_PHD_EXAMPLE_TRAJECTORY_GENERATION_HPP_

#include <array>
#include <vector>

struct GroundTruth {
  size_t time_steps;
  std::vector<double> orientation;
  std::vector<std::pair<double, double>> velocity;
  std::vector<std::pair<double, double>> size;
  std::vector<std::array<double, 4u>> rotation;
  std::vector<std::array<double, 2u>> center;
};

extern GroundTruth GetGroundTruth(void);

#endif  //  GM_PHD_EXAMPLE_TRAJECTORY_GENERATION_HPP_
