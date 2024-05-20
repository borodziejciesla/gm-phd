#ifndef INCLUD_ALIASES_HPP_
#define INCLUD_ALIASES_HPP_

#include <Eigen/Dense>

namespace mot {
template <size_t size>
using Vector = Eigen::Vector<float, size>;

template <size_t size>
using SquareMatrix = Eigen::Matrix<float, size, size>;

template <size_t size_row, size_t size_col>
using Matrix = Eigen::Matrix<float, size_row, size_col>;

using SensorPoseVector = Eigen::Vector<float, 3u>;

using SensorPoseMatrix = Eigen::Matrix<float, 3u, 3u>;
}  // namespace mot

#endif  //  INCLUD_ALIASES_HPP_
