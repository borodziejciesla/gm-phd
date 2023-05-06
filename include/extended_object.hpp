#ifndef GM_PHD_INCLUDE_EXTENDED_OBJECT_HPP_
#define GM_PHD_INCLUDE_EXTENDED_OBJECT_HPP_

#include "value_with_covariance.hpp"

namespace mot {
  struct ExtentState {
    float length;
    float width;
    float orientation;

    Eigen::Matrix<double, 3u, 3u> covariance = Eigen::Matrix<double, 3u, 3u>::Zero();
  }

  template <size_t kinematic_state_size>
  struct ExtendedObject {
    ValueWithCovariance<kinematic_state_size> kinematic_state;
    ExtentState extent_state;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_EXTENDED_OBJECT_HPP_
