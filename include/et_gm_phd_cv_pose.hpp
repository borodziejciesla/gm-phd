#ifndef ET_GM_PHD_INCLUDE_GM_PHD_CV_HPP_
#define ET_GM_PHD_INCLUDE_GM_PHD_CV_HPP_

#include "et_gm_phd.hpp"

namespace mot {
  class EtGmPhdCvPose : public EtGmPhd<4u, 2u> {
    public:
      explicit EtGmPhdCvPose(const GmPhdCalibrations<4u, 2u> & calibrations);
      virtual ~EtGmPhdCvPose(void) = default;

    protected:
      std::pair<Eigen::Vector<double, 4u>, Eigen::Matrix<double, 4u, 4u>> 
        PredictHypothesis(const Eigen::Vector<double, 4u> & state, const Eigen::Matrix<double, 4u, 4u> & covariance);
      void PrepareTransitionMatrix(void);
      void PrepareProcessNoiseMatrix(void);

      void PredictBirths(void);
  };
} //  namespace mot

#endif  //  ET_GM_PHD_INCLUDE_GM_PHD_CV_HPP_
