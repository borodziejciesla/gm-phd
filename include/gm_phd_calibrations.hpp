#ifndef GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_
#define GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_

#include <array>
#include <functional>

#include "aliases.hpp"
#include "hypothesis.hpp"

namespace mot {
template <size_t state_size, size_t measurement_size>
struct GmPhdCalibrations {
  std::array<float, state_size> process_noise_diagonal = {};  // Process noise covariance matrix

  std::function<void(Hypothesis<state_size, measurement_size>&, const float)> predict_hypothesis;
  std::function<void(Hypothesis<state_size, measurement_size>&)> predict_observation;

  float pd = 0.8f;  // Probability of detection
  float ps = 0.8f;  // Probability of survival
  float kappa = 1.0e-9f;

  float truncation_threshold = 0.1f;
  float merging_threshold = 3.0f;
};
}  //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_
