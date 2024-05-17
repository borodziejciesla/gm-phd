#ifndef GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
#define GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_

#include <Eigen/Dense>

namespace mot {
  template <size_t size>
  struct ValueWithCovariance {
    Eigen::Vector<float, size> value = Eigen::Vector<float, size>::Zero();
    Eigen::Matrix<float, size, size> covariance = Eigen::Matrix<float, size, size>::Zero();
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
