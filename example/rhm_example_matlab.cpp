#include <string>

#include "sensor_data_reader.hpp"

#include "rhm_gm_phd.hpp"

int main() {
  const auto log_path = std::string("/home/maciej/Downloads/sensors_data.csv");
  mot::example::SensorDataReader radar_data_reader(log_path);

  for (auto radar_index = 0u; radar_index < radar_data_reader.GetScansNumber(); radar_index++) {
    const auto [timestamp, radar_scan] = radar_data_reader.GetNextScan();
  }

  return EXIT_SUCCESS;
}