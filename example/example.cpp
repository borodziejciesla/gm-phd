#include <algorithm>
#include <array>
#include <filesystem>
#include <initializer_list>
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
using Filter = gm_phd_example::GmPhdCv;
using Measurements = std::vector<Filter::Measurement>;

struct FilterParameters {
  double pd = 1.0;
  double ps = 0.95;
  double kappa =
      8.0 / ((clutter_area_max - clutter_area_min) * (clutter_area_max - clutter_area_min));
  double truncation_threshold = 0.01;
  double merging_threshold = 4.0;
  double extraction_threshold = 0.5;
};

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

struct ExperimentResult {
  std::vector<std::vector<Filter::Object>> objects_by_scan;
  std::vector<int> object_count;
  std::vector<int> truth_object_count;
  std::vector<double> weights_sum;
  std::vector<double> ospa_scores;
  double mean_ospa = 0.0;
};

Simulation GenerateSimulation(const unsigned int seed) {
  std::mt19937 random_engine(seed);
  std::uniform_int_distribution<size_t> objects_number_dist(3u, 7u);
  std::uniform_int_distribution<size_t> start_dist(0u, measurements_number - 20u);
  std::uniform_int_distribution<size_t> duration_dist(20u, measurements_number);
  std::uniform_real_distribution<double> pose_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> velocity_dist(-1.0, 1.0);
  std::uniform_real_distribution<double> clutter_position_dist(clutter_area_min, clutter_area_max);
  std::poisson_distribution<size_t> clutter_count_dist(8.0);
  std::normal_distribution<double> measurement_noise(0.0, 0.3);
  const auto draw = [&random_engine](auto& distribution) { return distribution(random_engine); };
  Simulation simulation{std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number),
                        std::vector<Measurements>(measurements_number)};
  const auto objects_number = draw(objects_number_dist);
  std::vector<SimulatedObject> objects;
  objects.reserve(objects_number);

  for (size_t index = 0u; index < objects_number; ++index) {
    const auto first_scan = draw(start_dist);
    const auto last_scan = std::min(measurements_number - 1u, first_scan + draw(duration_dist));
    SimulatedObject object{first_scan, last_scan, Filter::StateVector::Zero()};
    object.state << draw(pose_dist), draw(velocity_dist), draw(pose_dist), draw(velocity_dist);
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
      measurement.value(0u) += draw(measurement_noise);
      measurement.value(1u) += draw(measurement_noise);
      measurement.covariance = 0.09 * Filter::MeasurementMatrix::Identity();
      simulation.target_measurements.at(scan).push_back(measurement);
      simulation.measurements.at(scan).push_back(measurement);
    }

    for (size_t clutter_index = 0u; clutter_index < draw(clutter_count_dist); ++clutter_index) {
      Filter::Measurement clutter;
      clutter.value << draw(clutter_position_dist), draw(clutter_position_dist);
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

ExperimentResult RunExperiment(const Simulation& simulation, const FilterParameters& parameters,
                               const bool print_scan_metrics = false) {
  Filter filter;
  const auto calibration = filter.SetCalibrations(
      "pd", parameters.pd, "ps", parameters.ps, "kappa", parameters.kappa, "truncation_threshold",
      parameters.truncation_threshold, "merging_threshold", parameters.merging_threshold,
      "extraction_threshold", parameters.extraction_threshold);
  if (!calibration) {
    throw std::runtime_error("Failed to set filter calibrations: " + calibration.error());
  }

  ExperimentResult result;
  result.objects_by_scan.reserve(measurements_number);
  result.object_count.reserve(measurements_number);
  result.truth_object_count.reserve(measurements_number);
  result.weights_sum.reserve(measurements_number);
  result.ospa_scores.reserve(measurements_number);
  for (auto index = 0u; index < measurements_number; index++) {
    filter.Run(static_cast<double>(index + 1u) * dt, simulation.measurements.at(index));
    const auto& objects = filter.GetObjects();
    result.objects_by_scan.push_back(objects);
    result.object_count.push_back(static_cast<int>(objects.size()));
    result.truth_object_count.push_back(static_cast<int>(simulation.truth.at(index).size()));
    result.weights_sum.push_back(filter.GetWeightsSum());
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
    result.ospa_scores.push_back(
        mot::OspaDistance(estimated_positions, true_positions, ospa_cutoff, ospa_order));
    if (print_scan_metrics) {
      std::cout << "Scan index: " << index << ", Objects number: " << objects.size()
                << ", Truth objects: " << simulation.truth.at(index).size()
                << ", Weights Sum: " << filter.GetWeightsSum()
                << ", OSPA: " << result.ospa_scores.back() << " m\n";
    }
  }
  result.mean_ospa = std::accumulate(result.ospa_scores.begin(), result.ospa_scores.end(), 0.0) /
                     result.ospa_scores.size();
  return result;
}

double EvaluateParameters(const FilterParameters& parameters,
                          const std::vector<Simulation>& simulations) {
  auto total_ospa = 0.0;
  for (const auto& simulation : simulations) {
    total_ospa += RunExperiment(simulation, parameters).mean_ospa;
  }
  return total_ospa / simulations.size();
}

FilterParameters TuneParameters(const std::vector<Simulation>& simulations) {
  FilterParameters best;
  auto best_score = EvaluateParameters(best, simulations);
  std::cout << "Initial mean OSPA: " << best_score << " m\n";

  auto tune_parameter = [&](const std::string& name, double FilterParameters::* parameter,
                            const std::initializer_list<double> candidates) {
    auto parameter_best = best;
    auto parameter_best_score = best_score;
    for (const auto value : candidates) {
      auto candidate = best;
      candidate.*parameter = value;
      const auto score = EvaluateParameters(candidate, simulations);
      std::cout << "  " << name << '=' << value << " -> mean OSPA " << score << " m\n";
      if (score < parameter_best_score) {
        parameter_best = candidate;
        parameter_best_score = score;
      }
    }
    best = parameter_best;
    best_score = parameter_best_score;
    std::cout << "  Selected " << name << '=' << best.*parameter << " (mean OSPA " << best_score
              << " m)\n";
  };

  tune_parameter("ps", &FilterParameters::ps, {0.8, 0.85, 0.90, 0.95, 0.99});
  tune_parameter("kappa", &FilterParameters::kappa, {0.005, 0.01, 0.02, 0.05});
  tune_parameter("truncation_threshold", &FilterParameters::truncation_threshold,
                 {0.005, 0.01, 0.02, 0.04, 0.08, 0.1});
  tune_parameter("merging_threshold", &FilterParameters::merging_threshold,
                 {0.5, 1.0, 1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 8.0});
  tune_parameter("extraction_threshold", &FilterParameters::extraction_threshold,
                 {0.3, 0.4, 0.5, 0.6, 0.7, 0.8});
  return best;
}

void SavePlots(const Simulation& simulation, const ExperimentResult& result) {
  const auto output_directory = std::filesystem::current_path() / "fig";
  std::filesystem::create_directories(output_directory);
  const auto estimated_trajectories = BuildEstimatedTrajectories(result.objects_by_scan);
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
  PlotSeries(scan_times,
             std::vector<double>(result.object_count.begin(), result.object_count.end()),
             "Estimated objects", "blue", ".");
  PlotSeries(
      scan_times,
      std::vector<double>(result.truth_object_count.begin(), result.truth_object_count.end()),
      "Ground truth", "red", "", "-");
  SaveAndShow(output_directory / "object_count.png");

  plt::figure_size(1200, 780);
  plt::title("PHD weight sum vs ground truth object count");
  plt::xlabel("Time [s]");
  plt::ylabel("Sum of hypothesis weights [-]");
  PlotSeries(
      scan_times,
      std::vector<double>(result.truth_object_count.begin(), result.truth_object_count.end()),
      "Ground truth object count", "red", "", "--");
  PlotSeries(scan_times, result.weights_sum, "PHD weight sum", "blue", "o", "-");
  SaveAndShow(output_directory / "phd_weight_sum.png");

  plt::figure_size(1200, 780);
  plt::title("OSPA distance (c=10 m, p=2)");
  plt::xlabel("Time [s]");
  plt::ylabel("OSPA [m]");
  PlotSeries(scan_times, result.ospa_scores, "OSPA", "blue", "", "-");
  SaveAndShow(output_directory / "ospa.png");
}

int main() {
  constexpr size_t tuning_runs = 50u;
  std::vector<Simulation> simulations;
  simulations.reserve(tuning_runs);
  for (size_t run = 0u; run < tuning_runs; ++run) {
    const auto seed = random_seed + static_cast<unsigned int>(run);
    std::cout << "Generated tuning run with seed: " << seed << '\n';
    simulations.push_back(GenerateSimulation(seed));
  }

  const auto best_parameters = TuneParameters(simulations);
  const auto best_result = RunExperiment(simulations.front(), best_parameters, true);
  std::cout << "Selected parameters: pd=" << best_parameters.pd << ", ps=" << best_parameters.ps
            << ", kappa=" << best_parameters.kappa
            << ", truncation_threshold=" << best_parameters.truncation_threshold
            << ", merging_threshold=" << best_parameters.merging_threshold
            << ", extraction_threshold=" << best_parameters.extraction_threshold << '\n'
            << "Mean OSPA on selected plot run (c=" << ospa_cutoff << " m, p=" << ospa_order
            << "): " << best_result.mean_ospa << " m\n"
            << "Mean OSPA across tuning runs: " << EvaluateParameters(best_parameters, simulations)
            << " m\n";
  SavePlots(simulations.front(), best_result);
  return 0;
}
