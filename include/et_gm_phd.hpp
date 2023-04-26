#ifndef GM_PHD_INCLUDE_ET_GM_PHD_HPP_
#define GM_PHD_INCLUDE_ET_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include <Eigen/Dense>

#include "gm_phd_calibrations.hpp"
#include "extended_object.hpp"

namespace mot {
  template <size_t state_size, size_t measurement_size>
  class EtGmPhd {
    public:
      using Object = ExtentState<state_size>;

    public:
      EtGmPhd() {}

      ~EtGmPhd(void) = default;

      const std::vector<Object> & GetObjects(void) const {
        return objects_;
      }

      void Run(const double timestamp, const std::vector<Measurement> & measurements) {
       //
      }

      double GetWeightsSum(void) const {
        return std::accumulate(hypothesis_.begin(), hypothesis_.end(),
          0.0,
          [](double sum, const Hypothesis & hypothesis) {
            return sum + hypothesis.weight;
          }
        );
      }

    protected:
      struct Hypothesis {
        Hypothesis(void) = default;
        Hypothesis(const Hypothesis&) = default;
        Hypothesis(Hypothesis&&) = default;
        Hypothesis & operator=(const Hypothesis&) = default;
        Hypothesis(const double w, const StateSizeVector s, const StateSizeMatrix c)
          : weight{w}
          , state{s}
          , covariance{c} {}

        bool operator==(const Hypothesis & arg) {
          return (weight == arg.weight)
            && (state == arg.state)
            && (covariance == arg.covariance);
        }

        double weight = 0.0;
        StateSizeVector state = StateSizeVector::Zero();
        StateSizeMatrix covariance = StateSizeMatrix::Zero();
      };
      
      std::vector<Object> objects_;

    private:
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_ET_GM_PHD_HPP_
