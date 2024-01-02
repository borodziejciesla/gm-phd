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

#include "distance_partitioner.hpp"
#include "gm_phd_calibrations.hpp"
#include "partition.hpp"
#include "value_with_covariance.hpp"

namespace mot {
  constexpr auto kinematic_state_size = 4u;
  constexpr auto extend_state_size = 3u;
  constexpr auto measurement_size = 2u;
  constexpr auto nx = extend_state_size + kinematic_state_size;
  constexpr auto augmented_size = nx + 1u + measurement_size;

  class RhmGmPhd {
    public:
      using KinematicStateSizeVector = Eigen::Vector<double, kinematic_state_size>;
      using KinematicStateSizeMatrix = Eigen::Matrix<double, kinematic_state_size, kinematic_state_size>;
      using ExtendStateSizeVector = Eigen::Vector<double, extend_state_size>;
      using ExtendStateSizeMatrix = Eigen::Matrix<double, extend_state_size, extend_state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;
      using SensorPoseVector = Eigen::Vector<double, 3u>;
      using SensorPoseMatrix = Eigen::Matrix<double, 3u, 3u>;
      using AugmentedVector = Eigen::Vector<double, augmented_size>;
      using AugmentedMatrix = Eigen::Matrix<double, augmented_size, augmented_size>;

      using UkfVector = Eigen::Vector<double, extend_state_size + 3u>;
      using UkfMatrix = Eigen::Matrix<double, extend_state_size + 3u, extend_state_size + 3u>;

      using ObservationMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;
      using TransitionMatrix = Eigen::Matrix<double, kinematic_state_size, kinematic_state_size>;

      using Measurement = ValueWithCovariance<measurement_size>;

      using Hypothesis = struct {
        double weight = 0.0;
        
        ValueWithCovariance<kinematic_state_size> kinematic;
        ValueWithCovariance<extend_state_size> extend;

        double s = 0.0;
        double var_s = 0.0;
        double var_zero = 0.0;
      };
      using HypothesisList = std::vector<Hypothesis>;

      using PartitioningList = std::vector<Partition>;

      using Object = struct {
        ValueWithCovariance<kinematic_state_size> kinematic;
        ValueWithCovariance<extend_state_size> extend;
      };
      using ObjectsList = std::vector<Object>;

      using UtPoints = struct {
        std::vector<UkfVector> x_pts;
        std::vector<double> wm;
        std::vector<double> wc;
        size_t points_number = 0u;
        std::vector<double> z_pts;
      };

    public:
      RhmGmPhd(void);
      ~RhmGmPhd(void);

      void Run(const double timestamp, const std::vector<Measurement>& measurements);
      const ObjectsList& GetObjects(void) const;

    private:
      void MakeBirth(void);
      void MakePrediction(void);
      void MakePartitioning(const std::vector<Measurement>& measurements);
      void MakeCorrection(const std::vector<Measurement>& measurements);
      void MakeKinematicCorrection(const Measurement& measurement);
      void MakeExtendCorrection(const Hypothesis& hypothesis, const Measurement& measurement);
      void MakePruning(void);
      void MakeMerging(void);
      void ExtractObjects(void);
      UtPoints GetUkfPointsAndWeights(const UkfVector& mu_a, const UkfMatrix& c_a, const double theta, const Measurement& measurement);
      UkfVector UtCorrectionStep(const UtPoints& ut_points);

      double prev_timestamp_ = 0.0;
      double time_delta_ = 0.0;

      HypothesisList hypothesis_list_;
      HypothesisList predicted_hypothesis_list_;

      PartitioningList partitions_;
      
      ObjectsList objects_list_;

      TransitionMatrix transition_matrix_;
      ObservationMatrix observation_matrix_;

      KinematicStateSizeMatrix q_;

      Eigen::Vector<double, nx> mu_ = Eigen::Vector<double, nx>::Zero();
      Eigen::Matrix<double, nx, nx> px_ = Eigen::Matrix<double, nx, nx>::Zero();
      Eigen::Vector2d m_ = Eigen::Vector2d::Zero();
      double pl_ = 1.0;
      double scale_ = 0.7;
      double scale_var_ = 0.08;

      double gamma_ = 1e-5;
      double pd_ = 0.9;
      double ps_ = 0.9;

      const double alpha_ = 1.0;
      const double beta_ = 2.0;
      const double kappa_ = 0.0;

      double s_ = 0.0;
      double z_mean_ = 0.0;
      double z_mean_var_ = 0.0;

      int jk_ = 0u;

      Eigen::Matrix<double, 2u, 4u> h_ = Eigen::Matrix<double, 2u, 4u>::Zero();

      DistancePartitioner<Measurement> distance_partitioner_;
      std::vector<Hypothesis> pruned_hypothesis_;

      /* B = 0.5 * (A + A') */
      template <size_t matrix_size>
      Eigen::Matrix<double, matrix_size, matrix_size> MakeMatrixSymetric(Eigen::Matrix<double, matrix_size, matrix_size> & matrix) {
        Eigen::Matrix<double, matrix_size, matrix_size> output;
        
        output = 0.5 * (matrix + matrix.transpose());

        return output;
      }
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_RHM_GM_PHD_HPP_
