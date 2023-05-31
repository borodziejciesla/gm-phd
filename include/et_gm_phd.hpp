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
#include "partitioning.hpp"

namespace mot {
  template <size_t state_size, size_t measurement_size>
  class EtGmPhd {
    public:
      using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Measurement = ValueWithCovariance<measurement_size>;
      using Object = ExtendedObject<state_size>;

    public:
      EtGmPhd(const GmPhdCalibrations<4u, 2u> & calibrations)
        : calibrations_{calibrations} {}

      ~EtGmPhd(void) = default;

      const std::vector<Object> & GetObjects(void) const {
        return objects_;
      }

      void Run(const double timestamp, const std::vector<Measurement> & measurements) {
        SetTimestamps(timestamp);
        // Make partitioning - create object input hypothesis
        PrepareInputHypothesis(measurements);
        // Run Filter
        Predict();
        UpdateMeasurements(measurements);
        // Post Processing
        Prune();
        ExtractObjects();
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
        Hypothesis(const double w, const StateSizeVector s, const StateSizeMatrix c, const ExtentState e)
          : weight{w}
          , state{s}
          , covariance{c}
          , extent_state{e} {}

        bool operator==(const Hypothesis & arg) {
          return (weight == arg.weight)
            && (state == arg.state)
            && (covariance == arg.covariance);
        }

        double weight = 0.0;
        StateSizeVector state = StateSizeVector::Zero();
        StateSizeMatrix covariance = StateSizeMatrix::Zero();
        ExtentState extent_state = ExtentState();
      };

      struct InputHypothesis {
        InputHypothesis(void) = default;
        InputHypothesis(const InputHypothesis&) = default;
        InputHypothesis(InputHypothesis&&) = default;
        InputHypothesis & operator=(const InputHypothesis&) = default;
        InputHypothesis(const double w, const StateSizeVector s, const StateSizeMatrix c, const ExtentState e)
          : weight{w}
          , state{s}
          , covariance{c} {}

        bool operator==(const Hypothesis & arg) {
          return (weight == arg.weight)
            && (state == arg.state)
            && (covariance == arg.covariance);
        }

        void AddDetectionIndex(const uint32_t detection_index) {
          associated_measurements_indices.push_back(detection_index);
        }

        double weight = 0.0;
        StateSizeVector state = StateSizeVector::Zero();
        StateSizeMatrix covariance = StateSizeMatrix::Zero();
        std::vector<size_t> associated_measurements_indices;
      };

      virtual Hypothesis PredictHypothesis(const Hypothesis & object) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;
      virtual void PredictBirths(void) = 0;

      double time_delta = 0.0;
      GmPhdCalibrations<state_size, measurement_size> calibrations_;
      std::vector<Hypothesis> predicted_hypothesis_;
      StateSizeMatrix transition_matrix_ = StateSizeMatrix::Zero();
      StateSizeMatrix process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

    private:
      void SetTimestamps(const double timestamp) {
        if (prev_timestamp_ != 0.0)
          time_delta = timestamp - prev_timestamp_;
        prev_timestamp_ = timestamp;
      }

      void Predict(void) {
        if (is_initialized_)
          predicted_hypothesis_.clear();
        else
          is_initialized_ = true;
        // PredictBirths();
        PredictExistingTargets();
      }

      void PredictExistingTargets(void) {
        // Prepare for prediction 
        PrepareTransitionMatrix();
        PrepareProcessNoiseMatrix();
        // Predict
        std::transform(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(predicted_hypothesis_),
          [this](const Hypothesis & hypothesis) {
            // Predict kinematics
            auto predicted_state = PredictHypothesis(hypothesis);
            // Predict weight
            predicted_state.weight *= calibrations_.ps;

            // Predict shape

            return Hypothesis(predicted_state);
          }
        );
      }

      void UpdateMeasurements(const std::vector<Measurement> & measurement) {
        // Clear hypothesis list
        hypothesis_.clear();
        // Create new hypothesis
        UpdateExistingHypothesis();
        MakeMeasurementUpdate(measurement);
      }

      void UpdateExistingHypothesis(void) {
        std::transform(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          std::back_inserter(hypothesis_),
          [this](const Hypothesis & hypothesis) {
            return Hypothesis((1.0 - (1.0 - std::exp(-gamma_)) * calibrations_.pd) * hypothesis.weight,
              hypothesis.state,
              hypothesis.covariance,
              hypothesis.extent_state);
          }
        );
      }

      void MakeMeasurementUpdate(const std::vector<Measurement> & measurements) {
        for (const auto & input_hypothesis : input_hypothesis_) {
          for (const auto & predicted_hypothesis : predicted_hypothesis_) {
            // Calculate Weight
            const auto wp = CalculateOmegaP(input_hypothesis);
            const auto gamma = CalculateGamma(input_hypothesis);
            const auto dw = CalculateDw(input_hypothesis, measurements);
            const auto phi_w = CalculateInputHypothesisPhi(predicted_hypothesis, input_hypothesis, measurements);

            const auto new_weight = wp * (gamma * calibrations_.pd / dw) * phi_w * (predicted_hypothesis.weight / calibrations_.ps);

            // Calculate kinematic state and covariance
            const auto [state, covariance] = UpdateKinematic(input_hypothesis, predicted_hypothesis, measurements);

            // Add element
            Hypothesis h;
            h.weight = new_weight;
            h.state = state;
            h.covariance = covariance;

            hypothesis_.push_back(h);
          }
        }
      }

      std::pair<StateSizeVector, StateSizeMatrix> UpdateKinematic(const InputHypothesis & input_hypothesis, const Hypothesis & predicted_hypothesis, const std::vector<Measurement> & measurements) {
        const auto kalman_gain = CalculateKalmaGain(input_hypothesis, predicted_hypothesis, measurements);
        const auto [state, covariance] = CalculateUpdatedKinematic(kalman_gain, input_hypothesis, predicted_hypothesis, measurements);

        return std::make_pair(state, covariance);
      }

      Eigen::MatrixXd CalculateKalmaGain(const InputHypothesis & input_hypothesis, const Hypothesis & predicted_hypothesis, const std::vector<Measurement> & measurements) const {
        const auto detections_number = input_hypothesis.associated_measurements_indices.size();
        
        Eigen::MatrixXd innovation_matrix(detections_number * measurement_size, detections_number * measurement_size);
        Eigen::MatrixXd observation_matrix(state_size, detections_number * measurement_size);
        
        for (auto row_index = 0u; row_index < detections_number; row_index++) {
          observation_matrix.block(0u, row_index * measurement_size, state_size, measurement_size) = calibrations_.observation_matrix.transpose();
          for (auto col_index = 0u; col_index < detections_number; col_index++) {
            const auto hph = calibrations_.observation_matrix * predicted_hypothesis.covariance * calibrations_.observation_matrix.transpose();
            if (row_index != col_index) {
              innovation_matrix.block(row_index * measurement_size, col_index * measurement_size, measurement_size, measurement_size) = hph;
            } else {
              innovation_matrix.block(row_index * measurement_size, col_index * measurement_size, measurement_size, measurement_size)
                = hph + measurements.at(row_index).covariance;
            }
          }
        }
        const auto innovation_matrix_inversed = innovation_matrix.inverse();
        
        return predicted_hypothesis.covariance * observation_matrix * innovation_matrix_inversed;
      }

      std::pair<StateSizeVector, StateSizeMatrix> CalculateUpdatedKinematic(const Eigen::MatrixXd & kalman_gain, const InputHypothesis & input_hypothesis, const Hypothesis & predicted_hypothesis, const std::vector<Measurement> & measurements) {
        const auto detections_number = input_hypothesis.associated_measurements_indices.size();

        Eigen::MatrixXd observation_matrix(detections_number * measurement_size, state_size);
        Eigen::MatrixXd observation(detections_number * measurement_size, 1u);

        for (auto index = 0u; index < detections_number; index++) {
          observation_matrix.block(index * measurement_size, 0u, measurement_size, state_size) = calibrations_.observation_matrix;
          observation.block(index * measurement_size, 0u, measurement_size, 1u) = measurements.at(input_hypothesis.associated_measurements_indices.at(index)).value;
        }

        const auto state = predicted_hypothesis.state + kalman_gain * (observation - observation_matrix * predicted_hypothesis.state);
        const auto covariance = (StateSizeMatrix::Identity() - kalman_gain * observation_matrix) * predicted_hypothesis.covariance;

        return std::make_pair(state, covariance);
      }

      double CalculateGamma(const InputHypothesis & input_hypothesis) const {
        return std::exp(-gamma_) * std::pow(gamma_, input_hypothesis.associated_measurements_indices.size());
      }

      double CalculateMeasurementPhi(const Measurement & measurement, const Hypothesis & hypothesis) const {
        return NormPdf(measurement.value,
          calibrations_.observation_matrix * hypothesis.state,
          measurement.covariance + calibrations_.observation_matrix * hypothesis.covariance * calibrations_.observation_matrix.transpose());
      }

      double CalculateInputHypothesisPhi(const Hypothesis & predicted_hypothesis, const InputHypothesis & input_hypothesis, const std::vector<Measurement> & measurements) const {
        double phi = 1.0;
        for (const auto det_index : input_hypothesis.associated_measurements_indices)
          phi *= CalculateMeasurementPhi(measurements.at(det_index), predicted_hypothesis);

        return phi / (lambda_ * ck_);
      }

      double CalculateOmegaP(const InputHypothesis & input_hypothesis) const {
        return 1.0; // Should be implemented in case of multiple partitioning
      }

      double CalculateDw(const InputHypothesis & input_hypothesis, const std::vector<Measurement> & measurements) const {
        auto dw = std::accumulate(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          0.0,
          [measurements, input_hypothesis, this](double sum, const Hypothesis & hypothesis) {
            const auto gamma = CalculateGamma(input_hypothesis);
            const auto phi_w = CalculateInputHypothesisPhi(hypothesis, input_hypothesis, measurements);
            return sum + (gamma * calibrations_.pd * phi_w * hypothesis.weight);
          }
        );

        if (predicted_hypothesis_.size() == 1u)
          dw += 1.0;
        return dw;
      }

      void Prune(void) {
        // Select elements with weigths over turncation threshold
        std::vector<Hypothesis> pruned_hypothesis;
        std::copy_if(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(pruned_hypothesis),
          [this](const Hypothesis & hypothesis) {
            return hypothesis.weight >= calibrations_.truncation_threshold;
          }
        );
        std::vector<std::pair<Hypothesis, bool>> pruned_hypothesis_marked;
        std::transform(pruned_hypothesis.begin(), pruned_hypothesis.end(),
          std::back_inserter(pruned_hypothesis_marked),
          [](const Hypothesis & hypothesis) {
            return std::make_pair(hypothesis, false);
          }        
        );

        // Merge hypothesis
        std::vector<Hypothesis> merged_hypothesis;
        auto non_marked_hypothesis_counter = [](size_t sum, const std::pair<Hypothesis, bool> & markable_hypothesis) {
          return sum + (markable_hypothesis.second ? 0u : 1u);
        };
        auto non_merged_hypothesis_number = std::accumulate(pruned_hypothesis_marked.begin(), pruned_hypothesis_marked.end(), 0u, non_marked_hypothesis_counter);

        while (non_merged_hypothesis_number > 0u) {
          auto I = pruned_hypothesis_marked | std::views::filter([](const std::pair<Hypothesis, bool> & hypothesis_mark) { return !hypothesis_mark.second; });

          // Select maximum weight element
          const auto maximum_weight_hypothesis = *std::max_element(I.begin(), I.end(),
            [](const std::pair<Hypothesis, bool> & a, const std::pair<Hypothesis, bool> & b) {
              return a.first.weight < b.first.weight;
            }
          );

          // Select hypothesis in merging threshold
          auto L = pruned_hypothesis_marked | std::views::filter(
            [maximum_weight_hypothesis,this](const std::pair<Hypothesis, bool> & markable_hypothesis) {
              const auto diff = markable_hypothesis.first.state - maximum_weight_hypothesis.first.state;
              const auto distance_matrix = diff.transpose() * markable_hypothesis.first.covariance.inverse() * diff;
              return (distance_matrix(0) < calibrations_.merging_threshold) && !markable_hypothesis.second;
            }
          );

          // Calculate new merged element
          const auto merged_weight = std::accumulate(L.begin(), L.end(),
            0.0,
            [](double sum, const std::pair<Hypothesis, bool> & hypothesis) {
              return sum + hypothesis.first.weight;
            }
          );
          
          StateSizeVector merged_state = StateSizeVector::Zero();
          for (const auto l : L)
            merged_state += (l.first.weight * l.first.state) / merged_weight;

          StateSizeMatrix merged_covariance = StateSizeMatrix::Zero();
          for (const auto l : L) {
            const auto diff = merged_state - l.first.state;
            merged_covariance += (l.first.covariance + diff * diff.transpose()) / merged_weight;
          }

          merged_hypothesis.push_back(Hypothesis(merged_weight, merged_state, merged_covariance, ExtentState()));
          // Remove L from I
          std::transform(L.begin(), L.end(),
           L.begin(),
            [](std::pair<Hypothesis, bool> & markable_hypothesis) {
              markable_hypothesis.second = true;
              return markable_hypothesis;
            }
          );
          //
          non_merged_hypothesis_number = std::accumulate(pruned_hypothesis_marked.begin(), pruned_hypothesis_marked.end(), 0u, non_marked_hypothesis_counter);
        }
        // Set final hypothesis
        hypothesis_ = merged_hypothesis;
      }

      void ExtractObjects(void) {
        objects_.clear();
        for (const auto & hypothesis : hypothesis_) {
          if (hypothesis.weight > 0.5) {
            Object object;
            object.kinematic_state.value = hypothesis.state;
            object.kinematic_state.covariance = hypothesis.covariance;
            objects_.push_back(object);
          }
        }
      }

      void PrepareInputHypothesis(const std::vector<Measurement> & measurements) {
        // Make partitioning
        const auto [cell_id, cell_numbers] = partitioner_.MakePartitioning(measurements);

        // Convert to input hypothesis
        input_hypothesis_.resize(cell_id);
        for (auto & ih : input_hypothesis_)
          ih.associated_measurements_indices.clear();

        for (auto detection_index = 0u; detection_index < cell_numbers.size(); detection_index++)
          input_hypothesis_.at(cell_numbers.at(detection_index)).associated_measurements_indices.push_back(detection_index);
      }

      static double NormPdf(const MeasurementSizeVector & z, const MeasurementSizeVector & nu, const MeasurementSizeMatrix & cov) {
        const auto diff = z - nu;
        const auto c = 1.0 / (std::sqrt(std::pow(std::numbers::pi, measurement_size) * cov.determinant()));
        const auto e = std::exp(-0.5 * diff.transpose() * cov.inverse() * diff);
        return c * e;
      }

      double prev_timestamp_ = 0.0;
      bool is_initialized_ = false;
      std::vector<Object> objects_;
      std::vector<InputHypothesis> input_hypothesis_;
      std::vector<Hypothesis> hypothesis_;
      DistancePartitioner<measurement_size> partitioner_;

      const double gamma_ = 1.0;
      const double lambda_ = 0.1;
      const double ck_ = 1.0e-9;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_ET_GM_PHD_HPP_
