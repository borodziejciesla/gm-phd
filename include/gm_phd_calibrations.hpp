#ifndef GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_
#define GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_

#include <array>

#include <Eigen/Dense>

namespace mot {
  template <size_t state_size, size_t measurement_size>
  struct GmPhdCalibrations {
    std::array<double, state_size> process_noise_diagonal = {};                 // Process noise covariance matrix

    Eigen::Matrix<double, measurement_size, state_size> observation_matrix;     // Observation matrix
    Eigen::Matrix<double, measurement_size, measurement_size> measurement_covariance;

    double pd = 0.95;                                                       // Probability of detection
    double ps = 0.8;                                                        // Probability of survival
    double kappa = 1.0e-9;

    double truncation_threshold = 0.01;
    double merging_threshold = 1.0;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_CALIBRATIONS_HPP_
