#ifndef GM_PHD_SRC_HELPERS_HPP_
#define GM_PHD_SRC_HELPERS_HPP_

#include <cmath>
#include <numbers>
#include <stdexcept>

#include <Eigen/Cholesky>
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
  if (!z.allFinite() || !nu.allFinite() || !cov.allFinite()) {
    throw std::domain_error("NormPdf requires finite inputs");
  }

  const auto diff = z - nu;
  const Eigen::LLT<Matrix<size>> factorization(cov);
  if (factorization.info() != Eigen::Success) {
    throw std::domain_error("NormPdf covariance must be positive definite (determinant " +
                            std::to_string(cov.determinant()) + ")");
  }

  const auto log_determinant =
      2.0 * factorization.matrixL().toDenseMatrix().diagonal().array().log().sum();
  const auto normalized_exponent =
      -0.5 * (static_cast<double>(size) * std::log(2.0 * std::numbers::pi) + log_determinant);
  const auto mahalanobis_distance = diff.dot(factorization.solve(diff));
  return std::exp(normalized_exponent - 0.5 * mahalanobis_distance);
}
}  //  namespace mot

#endif  //  GM_PHD_SRC_HELPERS_HPP_