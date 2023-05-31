#ifndef GM_PHD_INCLUDE_PARTITIONING_HPP_
#define GM_PHD_INCLUDE_PARTITIONING_HPP_

#include <vector>

#include "value_with_covariance.hpp"

namespace mot {
  template <uint8_t measurement_size>
  class DistancePartitioner {
    private:
      static constexpr int32_t invalid_index = -1;

    public:
      using Measurement = ValueWithCovariance<measurement_size>;

    public:
      DistancePartitioner(void) = default;
      virtual ~DistancePartitioner(void) = default;

      std::pair<int32_t, std::vector<int32_t>> MakePartitioning(const std::vector<Measurement> & measurements) {
        // Prepare distance matrix
        CalculateDistances(measurements);
        
        // Prepare cells number
        cell_numbers_.resize(measurements.size());
        std::fill(cell_numbers_.begin(), cell_numbers_.end(), invalid_index);

        // Main partitioning loop
        RunMainPartitioning(measurements);

        return std::pair<int32_t, std::vector<int32_t>>(cell_id_, cell_numbers_);
      }

    private:
      using DistanceMatrix = std::vector<std::vector<float>>;

    private:
      void CalculateDistances(const std::vector<Measurement> & measurements) {
        // Clear matrix
        for (auto & row : distance_matrix_)
          row.clear();
        distance_matrix_.clear();

        // Allocate new matrix
        distance_matrix_ = DistanceMatrix(measurements.size());
        for (auto & row : distance_matrix_)
          row.resize(measurements.size());
        
        // Calculate distances
        for (auto row_index = 0u; row_index < measurements.size(); row_index++) {
          for (auto col_index = row_index; col_index < measurements.size(); col_index++) {
            const auto distance = CalculateMahalanobisDistance(measurements.at(row_index), measurements.at(col_index));
            distance_matrix_.at(row_index).at(col_index) = distance;
            distance_matrix_.at(col_index).at(row_index) = distance;
          }
        }
      }

      static float CalculateMahalanobisDistance(const Measurement & z_i, const Measurement & z_j) {
        const auto diff = z_i.value - z_j.value;
        const auto covariance = z_i.covariance + z_j.covariance;

        const auto distance_raw = diff.transpose() * covariance.inverse() * diff;
        return distance_raw(0u);
      }

      void FindNeihgbours(const uint32_t i, const std::vector<Measurement> & measurements, const uint32_t cell_id) {
        for (auto j = 0u; j < measurements.size(); j++) {
          const auto is_different_index = (j != i);
          const auto is_in_maximum_range = (distance_matrix_.at(i).at(j) <= 100.0);
          const auto is_non_initialized = (cell_numbers_.at(j) == invalid_index);

          if (is_different_index && is_in_maximum_range && is_non_initialized) {
            cell_numbers_.at(j) = cell_id;
            FindNeihgbours(j, measurements, cell_id);
          }
        }
      }

      void RunMainPartitioning(const std::vector<Measurement> & measurements) {
        cell_id_ = 0u;
        for (auto i = 0; i < measurements.size(); i++) {
          if (cell_numbers_.at(i) == invalid_index) {
            cell_numbers_.at(i) = cell_id_;
            FindNeihgbours(i, measurements, cell_id_);
            cell_id_++;
          }
        }
      }

      DistanceMatrix distance_matrix_;
      std::vector<int32_t> cell_numbers_;
      int32_t cell_id_ = 0;
  };
} // namespace mot

#endif  //  GM_PHD_INCLUDE_PARTITIONING_HPP_
