#ifndef GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
#define GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_

#include "aliases.hpp"

namespace mot {
template <size_t size>
struct ValueWithCovariance {
  Vector<size> value = Vector<size>::Zero();
  SquareMatrix<size> covariance = SquareMatrix<size>::Zero();
};
}  //  namespace mot

#endif  //  GM_PHD_INCLUDE_VALUE_WITH_COVARIANCE_HPP_
