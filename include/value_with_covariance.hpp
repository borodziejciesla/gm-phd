#ifndef GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
#define GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_

#include <Eigen/Dense>

namespace mot {
  template <size_t size>
  struct ValueWithCovariance {
    Eigen::Vector<double, size> value = Eigen::Vector<double, size>::Zero();
    Eigen::Matrix<double, size, size> covariance = Eigen::Matrix<double, size, size>::Zero();
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
