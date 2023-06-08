#ifndef GM_PHD_INCLUDE_GM_PHD_CV_HPP_
#define GM_PHD_INCLUDE_GM_PHD_CV_HPP_

#include "gm_phd.hpp"

namespace mot {
  class GmPhdCvPose : public GmPhd<4u, 2u> {
    public:
      explicit GmPhdCvPose(const GmPhdCalibrations<4u, 2u> & calibrations);
      virtual ~GmPhdCvPose(void) = default;

    protected:
      std::pair<Eigen::Vector<double, 4u>, Eigen::Matrix<double, 4u, 4u>> 
        PredictHypothesis(const Eigen::Vector<double, 4u> & state, const Eigen::Matrix<double, 4u, 4u> & covariance);
      void PrepareTransitionMatrix(void);
      void PrepareProcessNoiseMatrix(void);

      void PredictBirths(void);
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_GM_PHD_CV_HPP_
