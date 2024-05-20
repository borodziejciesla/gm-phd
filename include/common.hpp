#ifndef INCLUDE_COMMON_HPP_
#define INCLUDE_COMMON_HPP_

#include <cmath>
#include <numeric>

#include "aliases.hpp"

namespace mot {
template <size_t size>
float NormPdf(const Vector<size>& z, const Vector<size>& nu, const SquareMatrix<size>& cov) {
  const auto diff = z - nu;
  const auto c = 1.0f / (std::sqrt(std::pow(std::numbers::pi, size) * cov.determinant()));
  const auto e = std::exp(-0.5f * diff.transpose() * cov.inverse() * diff);
  return c * e;
}
}  //  namespace mot

#endif  //  INCLUDE_COMMON_HPP_
