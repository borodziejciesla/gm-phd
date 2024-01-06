#include "sensor_data_reader.hpp"

#include <cmath>

namespace mot::example {
  SensorDataReader::SensorDataReader(const std::string& file_path) {
    reader_ = std::make_unique<csv::CSVReader>(file_path);

    RadarScan radar_scan;
    MeasurementStruct measurement;

    for (csv::CSVRow& row : *reader_) {
      const auto scan_index = row["scan_index"].get<size_t>();

      if ((previous_scan_index_ != scan_index) && (previous_scan_index_ != 0u)) {
        radar_scans_.push_back(radar_scan);
        radar_scan.second.clear();
        previous_scan_index_++;
      }

      // Read line
      const auto timestamp = row["timestamp"].get<double>();
      const auto range = row["range"].get<double>();
      const auto range_std = row["range_std"].get<double>();
      const auto azimuth = row["azimuth"].get<double>();
      const auto azimuth_std = row["azimuth_std"].get<double>();
      const auto range_rate = row["range_rate"].get<double>();
      const auto range_rate_std = row["range_rate_std"].get<double>();

      measurement.value(0u) = range * std::cos(azimuth);
      measurement.value(1u) = range * std::sin(azimuth);

      measurement.covariance(0u, 0u) = range * std::cos(azimuth);
      measurement.covariance(0u, 1u) = measurement.covariance(1u, 0u) = range * std::cos(azimuth);
      measurement.covariance(1u, 1u) = range * std::cos(azimuth);

      radar_scan.first = timestamp;
      radar_scan.second.push_back(measurement);
    }

    previous_scan_index_ = 0u;
  }

  RadarScan SensorDataReader::GetNextScan(void) {
    return radar_scans_.at(previous_scan_index_++);
  }
}