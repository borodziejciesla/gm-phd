#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "sensor_data_reader.hpp"

#include "rhm_gm_phd.hpp"

static std::vector<mot::example::MeasurementStruct> ConvertScsToVcs(const std::vector<mot::example::MeasurementStruct> & scs_measurements);

int main() {
  /************************** Create Reader Object ***************************/
  const auto log_path = std::string("/home/maciej/Downloads/sensors_data.csv");
  mot::example::SensorDataReader radar_data_reader(log_path);

  /************************** Define tracker object **************************/
  mot::RhmGmPhd::RhmGmPhdCalibrations calibrations;
  calibrations.process_noise_diagonal = {0.5, 0.5, 2.0, 2.0};
  calibrations.observation_matrix = Eigen::Matrix<double, 2u, 4u>::Zero();
  calibrations.observation_matrix(0u, 0u) = 1.0;
  calibrations.observation_matrix(1u, 1u) = 1.0;
  calibrations.measurement_covariance = 0.2 * Eigen::Matrix<double, 2u, 2u>::Identity();

  mot::RhmGmPhd gm_phd_filter = mot::RhmGmPhd(calibrations);

  // Run simulation
  for (auto radar_index = 0u; radar_index < radar_data_reader.GetScansNumber(); radar_index++) {
    const auto [timestamp, radar_scan_scs] = radar_data_reader.GetNextScan();
    const auto radar_scan_vcs = ConvertScsToVcs(radar_scan_scs);

    gm_phd_filter.Run(timestamp, radar_scan_vcs);

    const auto objects = gm_phd_filter.GetObjects();
    std::cout << "Objects number: " << objects.size() << std::endl;
  }

  return EXIT_SUCCESS;
}

/***********************************************************************/
std::vector<mot::example::MeasurementStruct> ConvertScsToVcs(const std::vector<mot::example::MeasurementStruct> & scs_measurements) {
  std::vector<mot::example::MeasurementStruct> vcs_measurements;
  std::transform(scs_measurements.begin(), scs_measurements.end(),
    std::back_inserter(vcs_measurements),
    [](const mot::example::MeasurementStruct& scs_detection) -> mot::example::MeasurementStruct {
      mot::example::MeasurementStruct vcs_detection;

      vcs_detection.value(0u) = scs_detection.value(0) * std::cos(scs_detection.value(1));
      vcs_detection.value(1u) = scs_detection.value(0) * std::sin(scs_detection.value(1));

      vcs_detection.covariance = 0.2 * Eigen::Matrix2d::Identity();

      return vcs_detection;
    }
  );

  return vcs_measurements;
}
