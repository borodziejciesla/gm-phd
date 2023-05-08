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
        MakeDistancePartitioning(measurements);
        // Run Filter
        Predict();
        Update(measurements);
        // Post Processing
        //Prune();
        //ExtractObjects();
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
            const auto predicted_state = PredictHypothesis(hypothesis);
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
            return Hypothesis((1.0 - (1.0 - std::exp(-gamma_))) * hypothesis.weight,
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
            const auto dw = CalculateDw();
            const auto phi_w = CalculateInputHypothesisPhi(input_hypothesis, measurements);

            const auto new_weight = wp * (gamma * calibrations_.pd / dw) * phi_w * (predicted_hypothesis.weight / calibrations.ps);

            // Calculate kinematic state and covariance
          }
        }
      }

      double CalculateGamma(const InputHypothesis & input_hypothesis) const {
        return ste::exp(-gamma_) * std::pow(gamma_, input_hypothesis.associated_measurements_indices.size());
      }

      double CalculateMeasurementPhi(const Measurement & measurement, const Hypothesis & hypothesis) const {
        return NormPdf(measurement.value,
          h_ * hypothesis.state,
          measurement.covariance + h_ * hypothesis.covariance * h_.transpose);
      }

      double CalculateInputHypothesisPhi(const InputHypothesis & input_hypothesis, const std::vector<Measurement> & measurements) const {
        double phi = 1.0;

        for (const auto det_index : input_hypothesis.associated_measurements_indices)
          phi *= (CalculateMeasurementPhi(input_hypothesis) / (measurements.at(det_index)));

        return phi;
      }

      double CalculateOmegaP(const InputHypothesis & input_hypothesis) const {
        return 1.0;
      }

      double CalculateDw(void) const {
        auto dw = std::accumulate(predicted_hypothesis_.begin(), predicted_hypothesis_.end(),
          0.0,
          [this](double sum, const Hypothesis & hypothesis) {
            return sum + 1.0;
          }
        );

        if (input_hypothesis.associated_measurements_indices.size() == 1u)
          dw += 1.0;

        return 1.0;
      }

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
            const auto distance = CalculateMahalanobisDistance(measurements.at(row_index), measurements.at(col_index));
            distance_matrix_.at(row_index).at(col_index) = distance;
            distance_matrix_.at(col_index).at(row_index) = distance;
          }
        }
      }

      void PrepareInputHypothesis(const std::vector<Measurement> & measurements) {
        input_hypothesis_.reshape(cell_id_);
        for (auto detection_index = 0u; detection_index < cell_numbers_.size(); detection_index++)
          input_hypothesis_.at(cell_numbers_.at(detection_index)) = detection_index;
      }

      void MakeDistancePartitioning(const std::vector<Measurement> & measurements) {
        // Prepare distance matrix
        CalculateDistances(measurements);
        
        // Prepare cells number
        cell_numbers_.resize(measurements.size());
        std::fill(cell_numbers_.begin(), cell_numbers_.end(), 0u);
        
        // Main partitioning loop
        cell_id_ = 0u;
        for (auto i = 0; i < measurements.size(); i++) {
          if (cell_numbers_.at(i) == 0u) {
            cell_numbers_.at(i) = cell_id_;
            FindNeihgbours(i, measurements, cell_id_);
            cell_id_++;
          }
        }
      }

      void FindNeihgbours(const uint32_t i, const std::vector<Measurement> & measurements, const uint32_t cell_id) {
        for (auto j = 0u; j < measurements.size(); j++) {
          const auto is_different_index = (j != i);
          const auto is_in_maximum_range = (distance_matrix_.at(i).at(j) <= 10.0);
          const auto is_non_initialized = (cell_numbers_.at(j) == 0u);

          if (is_different_index && is_in_maximum_range && is_non_initialized) {
            cell_numbers_.at(j) = cell_id;
            FindNeihgbours(j, measurements, cell_id);
          }
        }
      }

      static float CalculateMahalanobisDistance(const Measurement & z_i, const Measurement & z_j) {
        const auto diff = z_i.value - z_j.value;
        const auto covariance = z_i.covariance + z_j.covariance;

        const auto distance_raw = diff.transpose() * covariance.inverse() * diff;
        return distance_raw(0u);
      }

      static double NormPdf(const MeasurementSizeVector & z, const MeasurementSizeVector & nu, const MeasurementSizeMatrix & cov) {
        const auto diff = z - nu;
        const auto c = 1.0 / (std::sqrt(std::pow(std::numbers::pi, measurement_size) * cov.determinant()));
        const auto e = std::exp(-0.5 * diff.transpose() * cov.inverse() * diff);
        return c * e;
      }

      using DistanceMatrix = std::vector<std::vector<float>>;

      double prev_timestamp_ = 0.0;
      bool is_initialized_ = false;
      int32_t cell_id_ = 0u;
      DistanceMatrix distance_matrix_;
      std::vector<uint32_t> cell_numbers_;
      std::vector<Object> objects_;
      std::vector<InputHypothesis> input_hypothesis_;
      std::vector<Hypothesis> hypothesis_;

      const double gamma_ = 10.0;
      const double lambda = 1.0;
      const double ck = 1.0;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_ET_GM_PHD_HPP_
