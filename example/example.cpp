#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "cv_gm_phd.hpp"
#include "matplotlibcpp.hpp"
#include "ospa.hpp"

namespace plt = matplotlibcpp;

constexpr size_t measurements_number = 100u;
constexpr double dt = 0.1;
constexpr double clutter_area_min = -20.0;
constexpr double clutter_area_max = 20.0;
constexpr double ospa_cutoff = 10.0;
constexpr double ospa_order = 2.0;

constexpr unsigned int random_seed = 20261018u;
std::mt19937 e1(random_seed);

std::uniform_int_distribution<size_t> objects_number_dist(3u, 7u);
std::uniform_int_distribution<size_t> start_dist(0u, measurements_number - 20u);
std::uniform_int_distribution<size_t> duration_dist(20u, measurements_number);
std::uniform_real_distribution<double> pose_dist(-10.0, 10.0);
std::uniform_real_distribution<double> velocity_dist(-1.0, 1.0);
std::uniform_real_distribution<double> clutter_position_dist(clutter_area_min, clutter_area_max);
std::poisson_distribution<size_t> clutter_count_dist(8.0);
using Filter = gm_phd_example::GmPhdCv;
using Measurements = std::vector<Filter::Measurement>;

struct SimulatedObject {
  size_t first_scan;
  size_t last_scan;
  Filter::StateVector state;
};

struct Simulation {
  std::vector<Measurements> measurements;
  std::vector<Measurements> target_measurements;
  std::vector<Measurements> truth;
  std::vector<Measurements> clutter;
};

Simulation GenerateSimulation(void) {
  std::normal_distribution<double> measurement_noise(0.0, 0.3);
  Simulation simulation{std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number)};
  const auto objects_number = objects_number_dist(e1);
  std::vector<SimulatedObject> objects;
  objects.reserve(objects_number);

  for (size_t index = 0u; index < objects_number; ++index) {
    const auto first_scan = start_dist(e1);
    const auto last_scan = std::min(measurements_number - 1u, first_scan + duration_dist(e1));
    SimulatedObject object{first_scan, last_scan, Filter::StateVector::Zero()};
    object.state << pose_dist(e1), velocity_dist(e1), pose_dist(e1), velocity_dist(e1);
    objects.push_back(object);
  }

  for (size_t scan = 0u; scan < measurements_number; ++scan) {
    for (const auto& object : objects) {
      if (scan < object.first_scan || scan > object.last_scan) {
        continue;
      }
      const auto elapsed = static_cast<double>(scan - object.first_scan) * dt;
      Filter::Measurement truth;
      truth.value << object.state(0u) + elapsed * object.state(1u),
          object.state(2u) + elapsed * object.state(3u);
      simulation.truth.at(scan).push_back(truth);

      auto measurement = truth;
      measurement.value(0u) += measurement_noise(e1);
      measurement.value(1u) += measurement_noise(e1);
      measurement.covariance = 0.09 * Filter::MeasurementMatrix::Identity();
      simulation.target_measurements.at(scan).push_back(measurement);
      simulation.measurements.at(scan).push_back(measurement);
    }

    for (size_t clutter_index = 0u; clutter_index < clutter_count_dist(e1); ++clutter_index) {
      Filter::Measurement clutter;
      clutter.value << clutter_position_dist(e1), clutter_position_dist(e1);
      clutter.covariance = 0.09 * Filter::MeasurementMatrix::Identity();
      simulation.clutter.at(scan).push_back(clutter);
      simulation.measurements.at(scan).push_back(clutter);
    }
  }
  return simulation;
}

struct EstimatedTrajectory {
  std::vector<double> time;
  std::vector<double> x;
  std::vector<double> y;
  size_t last_scan = 0u;
};

std::vector<EstimatedTrajectory> BuildEstimatedTrajectories(
    const std::vector<std::vector<Filter::Object>>& objects_by_scan) {
  constexpr size_t max_missed_scans = 2u;
  constexpr double max_link_distance = 5.0;
  struct Candidate {
    double distance;
    size_t trajectory;
    size_t object;
  };

  std::vector<EstimatedTrajectory> trajectories;
  for (size_t scan = 0u; scan < objects_by_scan.size(); ++scan) {
    std::vector<Candidate> candidates;
    for (size_t trajectory = 0u; trajectory < trajectories.size(); ++trajectory) {
      const auto& previous = trajectories[trajectory];
      if (scan - previous.last_scan > max_missed_scans) {
        continue;
      }
      for (size_t object = 0u; object < objects_by_scan[scan].size(); ++object) {
        const auto& estimate = objects_by_scan[scan][object];
        const auto distance = std::hypot(estimate.value(0u) - previous.x.back(),
                                         estimate.value(2u) - previous.y.back());
        if (distance <= max_link_distance) {
          candidates.push_back({distance, trajectory, object});
        }
      }
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
      return lhs.distance < rhs.distance;
    });

    std::vector<bool> matched_trajectories(trajectories.size(), false);
    std::vector<bool> matched_objects(objects_by_scan[scan].size(), false);
    for (const auto& candidate : candidates) {
      if (matched_trajectories[candidate.trajectory] || matched_objects[candidate.object]) {
        continue;
      }
      auto& trajectory = trajectories[candidate.trajectory];
      const auto& estimate = objects_by_scan[scan][candidate.object];
      trajectory.time.push_back(static_cast<double>(scan + 1u) * dt);
      trajectory.x.push_back(estimate.value(0u));
      trajectory.y.push_back(estimate.value(2u));
      trajectory.last_scan = scan;
      matched_trajectories[candidate.trajectory] = true;
      matched_objects[candidate.object] = true;
    }

    for (size_t object = 0u; object < objects_by_scan[scan].size(); ++object) {
      if (matched_objects[object]) {
        continue;
      }
      const auto& estimate = objects_by_scan[scan][object];
      trajectories.push_back({{static_cast<double>(scan + 1u) * dt},
                              {estimate.value(0u)},
                              {estimate.value(2u)},
                              scan});
    }
  }
  return trajectories;
}

void AppendSeparatedTrajectories(const std::vector<EstimatedTrajectory>& trajectories,
                                 std::vector<double>& x, std::vector<double>& y,
                                 std::vector<double>* time = nullptr) {
  const auto separator = std::numeric_limits<double>::quiet_NaN();
  for (const auto& trajectory : trajectories) {
    x.insert(x.end(), trajectory.x.begin(), trajectory.x.end());
    y.insert(y.end(), trajectory.y.begin(), trajectory.y.end());
    if (time != nullptr) {
      time->insert(time->end(), trajectory.time.begin(), trajectory.time.end());
    }
    x.push_back(separator);
    y.push_back(separator);
    if (time != nullptr) {
      time->push_back(separator);
    }
  }
}

void PlotSeries(const std::vector<double>& x, const std::vector<double>& y,
                const std::string& label, const std::string& color, const std::string& marker = "",
                const std::string& line_style = "none") {
  std::map<std::string, std::string> style = {
      {"label", label}, {"color", color}, {"linestyle", line_style}};
  if (!marker.empty()) {
    style["marker"] = marker;
  }
  if (!plt::plot(x, y, style)) {
    throw std::runtime_error("Failed to add plot series: " + label);
  }
}

void SaveAndShow(const std::filesystem::path& output_path) {
  plt::legend();
  auto* pyplot = PyImport_ImportModule("matplotlib.pyplot");
  if (pyplot == nullptr) {
    PyErr_Print();
    throw std::runtime_error("Failed to load matplotlib.pyplot");
  }
  const auto path_string = output_path.string();
  auto* result = PyObject_CallMethod(pyplot, "savefig", "s", path_string.c_str());
  Py_DECREF(pyplot);
  if (result == nullptr) {
    PyErr_Print();
    throw std::runtime_error("Failed to save plot to " + path_string);
  }
  Py_DECREF(result);
  std::cout << "Saved plot: " << path_string << '\n';
  plt::show();
  plt::close();
}

int main() {
  Filter gm_phd_filter;
  constexpr double clutter_intensity =
      8.0 / ((clutter_area_max - clutter_area_min) * (clutter_area_max - clutter_area_min));
  gm_phd_filter.SetCalibrations("pd", 1.0, "ps", 0.95, "kappa", clutter_intensity,
                                "truncation_threshold", 0.01, "merging_threshold", 4.0,
                                "extraction_threshold", 0.5);
  const auto simulation = GenerateSimulation();
  std::cout << "Random seed: " << random_seed << '\n';

  // Run Filter
  std::vector<std::vector<Filter::Object>> objects_all;
  std::vector<int> objects_number;
  std::vector<int> truth_objects_number;
  std::vector<double> weights_sum;
  std::vector<double> ospa_scores;
  for (auto index = 0u; index < measurements_number; index++) {
    gm_phd_filter.Run(static_cast<double>(index + 1u) * dt, simulation.measurements.at(index));
    const auto objects = gm_phd_filter.GetObjects();
    objects_all.push_back(objects);
    objects_number.push_back(objects.size());
    truth_objects_number.push_back(simulation.truth.at(index).size());
    weights_sum.push_back(gm_phd_filter.GetWeightsSum());
    std::vector<mot::OspaPosition> estimated_positions;
    estimated_positions.reserve(objects.size());
    for (const auto& object : objects) {
      estimated_positions.push_back({object.value(0u), object.value(2u)});
    }
    std::vector<mot::OspaPosition> true_positions;
    true_positions.reserve(simulation.truth.at(index).size());
    for (const auto& target : simulation.truth.at(index)) {
      true_positions.push_back({target.value(0u), target.value(1u)});
    }
    ospa_scores.push_back(
        mot::OspaDistance(estimated_positions, true_positions, ospa_cutoff, ospa_order));
    std::cout << "Scan index: " << index << ", Objects number: " << objects.size()
              << ", Truth objects: " << simulation.truth.at(index).size()
              << ", Weights Sum: " << gm_phd_filter.GetWeightsSum()
              << ", OSPA: " << ospa_scores.back() << " m\n";
  }
  const auto mean_ospa =
      std::accumulate(ospa_scores.begin(), ospa_scores.end(), 0.0) / ospa_scores.size();
  std::cout << "Mean OSPA (c=" << ospa_cutoff << " m, p=" << ospa_order << "): " << mean_ospa
            << " m\n";

  const auto output_directory = std::filesystem::current_path() / "fig";
  std::filesystem::create_directories(output_directory);
  const auto estimated_trajectories = BuildEstimatedTrajectories(objects_all);
  std::vector<double> objects_x;
  std::vector<double> objects_y;
  std::vector<double> objects_time;
  AppendSeparatedTrajectories(estimated_trajectories, objects_x, objects_y, &objects_time);

  std::vector<double> ref_objects_x;
  std::vector<double> ref_objects_y;
  std::vector<double> ref_objects_time;
  for (auto index = 0u; index < measurements_number; index++) {
    for (const auto& trajectories : simulation.truth.at(index)) {
      ref_objects_x.push_back(trajectories.value(0u));
      ref_objects_y.push_back(trajectories.value(1u));
      ref_objects_time.push_back(static_cast<double>(index + 1u) * dt);
    }
  }

  std::vector<double> meas_objects_x;
  std::vector<double> meas_objects_y;
  std::vector<double> meas_objects_time;
  for (auto index = 0u; index < measurements_number; index++) {
    for (const auto& measurement : simulation.target_measurements.at(index)) {
      meas_objects_x.push_back(measurement.value(0u));
      meas_objects_y.push_back(measurement.value(1u));
      meas_objects_time.push_back(static_cast<double>(index + 1u) * dt);
    }
  }

  std::vector<double> clutter_x;
  std::vector<double> clutter_y;
  for (const auto& scan : simulation.clutter) {
    for (const auto& measurement : scan) {
      clutter_x.push_back(measurement.value(0u));
      clutter_y.push_back(measurement.value(1u));
    }
  }

  plt::figure_size(1200, 780);
  plt::title("Estimated objects and measurements");
  plt::xlabel("X [m]");
  plt::ylabel("Y [m]");
  PlotSeries(objects_x, objects_y, "Estimated trajectories", "blue", "+", "-");
  PlotSeries(ref_objects_x, ref_objects_y, "Ground truth", "red", ".");
  PlotSeries(meas_objects_x, meas_objects_y, "Target measurements", "green", ".");
  PlotSeries(clutter_x, clutter_y, "Clutter", "magenta", "x");
  SaveAndShow(output_directory / "positions_xy.png");

  // Plot on directions
  plt::figure_size(1200, 780);
  plt::title("X coordinate over time");
  plt::xlabel("Time [s]");
  plt::ylabel("X [m]");
  PlotSeries(objects_time, objects_x, "Estimated trajectories", "blue", "+", "-");
  PlotSeries(ref_objects_time, ref_objects_x, "Ground truth", "red", ".");
  PlotSeries(meas_objects_time, meas_objects_x, "Target measurements", "green", ".");
  SaveAndShow(output_directory / "x_vs_time.png");

  plt::figure_size(1200, 780);
  plt::title("Y coordinate over time");
  plt::xlabel("Time [s]");
  plt::ylabel("Y [m]");
  PlotSeries(objects_time, objects_y, "Estimated trajectories", "blue", "+", "-");
  PlotSeries(ref_objects_time, ref_objects_y, "Ground truth", "red", ".");
  PlotSeries(meas_objects_time, meas_objects_y, "Target measurements", "green", ".");
  SaveAndShow(output_directory / "y_vs_time.png");

  std::vector<double> scan_times;
  scan_times.reserve(measurements_number);
  for (auto index = 0u; index < measurements_number; ++index) {
    scan_times.push_back(static_cast<double>(index + 1u) * dt);
  }

  plt::figure_size(1200, 780);
  plt::title("Extracted count (blue) vs ground truth (red)");
  plt::xlabel("Time [s]");
  plt::ylabel("Extracted object count [-]");
  PlotSeries(scan_times, std::vector<double>(objects_number.begin(), objects_number.end()),
             "Estimated objects", "blue", ".");
  PlotSeries(scan_times,
             std::vector<double>(truth_objects_number.begin(), truth_objects_number.end()),
             "Ground truth", "red", "", "-");
  SaveAndShow(output_directory / "object_count.png");

  plt::figure_size(1200, 780);
  plt::title("PHD weight sum vs ground truth object count");
  plt::xlabel("Time [s]");
  plt::ylabel("Sum of hypothesis weights [-]");
  PlotSeries(scan_times,
             std::vector<double>(truth_objects_number.begin(), truth_objects_number.end()),
             "Ground truth object count", "red", "", "--");
  PlotSeries(scan_times, weights_sum, "PHD weight sum", "blue", "o", "-");
  SaveAndShow(output_directory / "phd_weight_sum.png");

  plt::figure_size(1200, 780);
  plt::title("OSPA distance (c=10 m, p=2)");
  plt::xlabel("Time [s]");
  plt::ylabel("OSPA [m]");
  PlotSeries(scan_times, ospa_scores, "OSPA", "blue", "", "-");
  SaveAndShow(output_directory / "ospa.png");

  return 0;
}
