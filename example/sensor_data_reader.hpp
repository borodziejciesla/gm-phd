#ifndef GM_PHD_EXAMPLE_SENSOR_DATA_READER_HPP_
#define GM_PHD_EXAMPLE_SENSOR_DATA_READER_HPP_

#include <memory>
#include <vector>
#include <string>

#include <csv.hpp>

#include "value_with_covariance.hpp"

namespace mot::example {
  using MeasurementStruct = ValueWithCovariance<2u>;
  using RadarScan = std::pair<double, std::vector<MeasurementStruct>>;

  class SensorDataReader {
    public:
      explicit SensorDataReader(const std::string& file_path);

      RadarScan GetNextScan(void);

    private:
      std::unique_ptr<csv::CSVReader> reader_;
      std::vector<RadarScan> radar_scans_;
      size_t previous_scan_index_ = 0u;
  };
} // namespace mot::example

#endif  //  GM_PHD_EXAMPLE_SENSOR_DATA_READER_HPP_
