#ifndef GM_PHD_SRC_HELPERS_HPP_
#define GM_PHD_SRC_HELPERS_HPP_

#include <Eigen/Dense>

namespace mot {
template <size_t size>
using Vector = Eigen::Vector<double, size>;

template <size_t size>
using Matrix = Eigen::Matrix<double, size, size>;

template <size_t state_size, size_t measurement_size>
using ObservationMatrix = Eigen::Matrix<double, measurement_size, state_size>;

template <size_t state_size, size_t measurement_size>
using KalmanGainMatrix = Eigen::Matrix<double, state_size, measurement_size>;

template <size_t size>
struct ValueWithCovariance {
  Vector<size> value = Vector<size>::Zero();
  Matrix<size> covariance = Matrix<size>::Zero();
};

using SensorPoseVector = Vector<3u>;
using SensorPoseMatrix = Matrix<3u>;

/*
 * Computes the normalized probability density function for a multivariate normal distribution.
 *
 * @param z The measurement vector.
 * @param nu The mean vector.
 * @param cov The covariance matrix.
 * @return The normalized PDF value.
 */
template <size_t size>
double NormPdf(const Vector<size>& z, const Vector<size>& nu, const Matrix<size>& cov) {
  const auto diff = z - nu;
  const auto c = 1.0 / (std::sqrt(std::pow(std::numbers::pi, size) * cov.determinant()));
  const auto e = std::exp(-0.5 * diff.transpose() * cov.inverse() * diff);
  return c * e;
}
}  //  namespace mot

#endif  //  GM_PHD_SRC_HELPERS_HPP_