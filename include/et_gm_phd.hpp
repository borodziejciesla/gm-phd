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
    using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Measurement = ValueWithCovariance<measurement_size>;
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
      void CalculateDistances(const std::vector<Measurement> & measurements) {
        // Clear matrix
        for (auto & row : distance_matrix_)
          row.clear();
        distance_matrix_.clear();

        // Allocate new matrix
        distance_matrix_ = DistanceMatrix(measurements.size());
        for (auto & row : distance_matrix_)
          row = std::vector(measurements.size());
        
        // Calculate distances
        for (auto row_index = 0u; row_index < measurements.size(); row_index++) {
          for (auto col_index = row_index; col_index < measurements.size(); col_index++) {
            const auto distance = CalculateMahalanobisDistance(measurements.at(row_indes), measurements.at(col_index));
            distance_matrix_.at(row_index).at(col_index) = distance;
            distance_matrix_.at(col_index).at(row_index) = distance;
          }
        }
      }

      void MakeDistancePartitioning(const std::vector<Measurement> & measurements) {
        // Prepare distance matrix
        CalculateDistances(measurements);
        // Main partitioning loop
      }

      void FindNeihgbours(const uint32_t i, const std::vector<Measurement> & measurements) {
        for (auto j = 0u; j <= measurements.size(); j++) {
          if ((j != i)) {
            //
          }
        }
      }

      static float CalculateMahalanobisDistance(const Measurement & z_i, const Measurement & z_j) {
        const auto diff = z_i.value - z_j.value;
        const auto covariance = z_i.covariance + z_j.covariance;

        const auto distance_raw = diff.transpose() * covariance.inverse() * diff;
        return distance_raw(0u);
      }

      using DistanceMatrix = std::vector<std::vector<float>>;

      DistanceMatrix distance_matrix_;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_ET_GM_PHD_HPP_
