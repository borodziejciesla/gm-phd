#ifndef GM_PHD_SRC_DISTANCE_PARTITIONER_HPP_
#define GM_PHD_SRC_DISTANCE_PARTITIONER_HPP_

#include <algorithm>
#include <vector>

#include "partition.hpp"

namespace mot {
  template <typename MeasurementType>
  class DistancePartitioner {
    public:
      DistancePartitioner(void) = default;
      ~DistancePartitioner(void) = default;

      std::vector<Partition> MakePartitioning(const std::vector<MeasurementType>& measurements) {
        CalculateDistances(measurements);
        // Line 1
        cell_number_.resize(measurements.size());
        std::fill(cell_number_.begin(), cell_number_.end(), 0u);
        // Line 2
        cell_id_ = 1u;
        // Line 3
        for (auto i = 0u; i < measurements.size(); i++) {
          // Line 4
          if (cell_number_.at(i) == 0u) {
            // Line 5
            cell_number_.at(i) = cell_id_;
            // Line 6
            FindNeighbours(i, measurements);
            // Line 7
            cell_id_++;
          }
        }

        // Convert
        ConvertToPartitions();
        return partitions_;
      }

    private:
      void CalculateDistances(const std::vector<MeasurementType>& measurements) {
        for (auto row_index = 0u; row_index < measurements.size(); row_index++) {
          for (auto col_index = 0u; col_index < measurements.size(); col_index++) {
            distances_[row_index][col_index] = (measurements[row_index].value - measurements[col_index].value).norm();
          }
        }
      }

      void FindNeighbours(const size_t i, const std::vector<MeasurementType>& measurements) {
        for (auto j = 0u; j < measurements.size(); j++) {
          if ((i != j) && (distances_[i][j] < threshold_) && (cell_number_.at(j) == 0u)) {
            cell_number_.at(j) = cell_id_;
            FindNeighbours(j, measurements);
          }
        }
      }

      void ConvertToPartitions(void) {
        partitions_.resize(cell_id_);

        for (auto index = 0u; index < cell_id_; index++) {
          partitions_.at(index).points.clear();
          for (auto i = 0u; i < cell_number_.size(); i++) {
            if (cell_number_.at(i) == index)
              partitions_.at(index).points.push_back(i);
          }
        }
      }

      std::vector<Partition> partitions_;
      std::vector<std::vector<double>> distances_;
      std::vector<size_t> cell_number_;
      const double threshold_ = 1.0;
      size_t cell_id_ = 0u;
  };
} //  namespace mot

#endif  //  GM_PHD_SRC_DISTANCE_PARTITIONER_HPP_
