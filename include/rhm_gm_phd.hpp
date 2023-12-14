#ifndef GM_PHD_INCLUDE_RHM_GM_PHD_HPP_
#define GM_PHD_INCLUDE_RHM_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "gm_phd_calibrations.hpp"
#include "value_with_covariance.hpp"

namespace mot {
  constexpr auto kinematic_state_size = 4u;
  constexpr auto extend_state_size = 3u;
  constexpr auto measurement_size = 2u;

  struct Partition {
    std::vector<int> points;
  };

  class RhmGmPhd {
    public:
      using KinematicStateSizeVector = Eigen::Vector<float, kinematic_state_size>;
      using KinematicStateSizeMatrix = Eigen::Matrix<float, kinematic_state_size, kinematic_state_size>;
      using ExtendStateSizeVector = Eigen::Vector<float, extend_state_size>;
      using ExtendStateSizeMatrix = Eigen::Matrix<float, extend_state_size, extend_state_size>;
      using MeasurementSizeVector = Eigen::Vector<float, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<float, measurement_size, measurement_size>;
      using SensorPoseVector = Eigen::Vector<float, 3u>;
      using SensorPoseMatrix = Eigen::Matrix<float, 3u, 3u>;

      using ObservationMatrix = Eigen::Matrix<float, measurement_size, measurement_size>;
      using TransitionMatrix = Eigen::Matrix<float, kinematic_state_size, kinematic_state_size>;

      using Measurement = ValueWithCovariance<measurement_size>;

      using Hypothesis = struct {
        float weight = 0.0f;
        
        ValueWithCovariance<kinematic_state_size> kinematic;
        ValueWithCovariance<extend_state_size> extend;
      };
      using HypothesisList = std::vector<Hypothesis>;

      using PartitioningList = std::vector<Partition>;

      using Object = struct {
        ValueWithCovariance<kinematic_state_size> kinematic;
        ValueWithCovariance<extend_state_size> extend;
      };
      using ObjectsList = std::vector<Object>;

    public:
      RhmGmPhd();
      ~RhmGmPhd(void);

      void Run(const double timestamp, const std::vector<Measurement>& measurements);
      const ObjectsList& GetObjects(void) const;

    private:
      void MakePrediction(void);
      void MakeCorrection(const std::vector<Measurement>& measurements);
      void MakePruning(void);
      void MakeMerging(void);
      void ExtractObjects(void);

      double prev_timestamp_ = 0.0;
      float time_delta_ = 0.0f;

      HypothesisList hypothesis_list_;
      HypothesisList predicted_hypothesis_list_;

      PartitioningList partitions_;
      
      ObjectsList objects_list_;

      TransitionMatrix transition_matrix_;
      ObservationMatrix observation_matrix_;

      KinematicStateSizeMatrix q_;

      float gamma_ = 1e-5f;
      float pd_ = 0.9f;
      float ps_ = 0.9f;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_RHM_GM_PHD_HPP_
