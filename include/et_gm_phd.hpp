#ifndef GM_PHD_INCLUDE_ET_GM_PHD_HPP_
#define GM_PHD_INCLUDE_ET_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include <Eigen/Dense>

#include "base_gm_phd.hpp"
#include "extended_object.hpp"
#include "gm_phd_calibrations.hpp"
#include "partitioning.hpp"

namespace mot {
  /* Hypothesis structure */
  typedef <size_t state_size>
  struct EtHypothesis {
    EtHypothesis(void) = default;
    EtHypothesis(const EtHypothesis&) = default;
    EtHypothesis(EtHypothesis&&) = default;
    EtHypothesis & operator=(const EtHypothesis&) = default;
    EtHypothesis(const double w, const Eigen::Vector<double, state_size> s, const Eigen::Matrix<double, state_size, state_size> c, const ExtentState e)
      : weight{w}
      , state{s}
      , covariance{c}
      , extent_state{e} {}

    bool operator==(const EtHypothesis & arg) {
      return (weight == arg.weight)
        && (state == arg.state)
        && (covariance == arg.covariance);
    }

    double weight = 0.0;
    Eigen::Vector<double, state_size> state = Eigen::Vector<double, state_size>::Zero();
    Eigen::Matrix<double, state_size, state_size> covariance = Eigen::Matrix<double, state_size, state_size>::Zero();
    ExtentState extent_state = ExtentState();
  };

  /* Predicted Hypothesis structure */
  typedef <size_t state_size, size_t measurement_size>
  struct EtInputHypothesis {
    EtInputHypothesis(void) = default;
    EtInputHypothesis(const EtInputHypothesis&) = default;
    EtInputHypothesis(EtInputHypothesis&&) = default;
    EtInputHypothesis & operator=(const EtInputHypothesis&) = default;
    EtInputHypothesis(const double w, const Eigen::Vector<double, state_size> s, const Eigen::Matrix<double, state_size, state_size> c, const ExtentState e)
      : weight{w}
      , state{s}
      , covariance{c} {}

    bool operator==(const EtHypothesis & arg) {
      return (weight == arg.weight)
        && (state == arg.state)
        && (covariance == arg.covariance);
    }

    void AddDetectionIndex(const uint32_t detection_index) {
      associated_measurements_indices.push_back(detection_index);
    }

    double weight = 0.0;
    Eigen::Vector<double, state_size> state = Eigen::Vector<double, state_size>::Zero();
    Eigen::Matrix<double, state_size, state_size> covariance = Eigen::Matrix<double, state_size, state_size>::Zero();
    std::vector<size_t> associated_measurements_indices;
  };

  /* ET-GM-PHD Class */
  template <size_t state_size, size_t measurement_size, typename Hypothesis = EtHypothesis<state_size>, typename PredictedHypothesis = EtInputHypothesis<state_size, measurement_size>>
  class EtGmPhd : public BaseGmPhd<state_size, measurement_size, EtHypothesis<state_size>, EtInputHypothesis<state_size, measurement_size>> {
    public:
      using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Measurement = ValueWithCovariance<measurement_size>;
      using Object = ExtendedObject<state_size>;

    public:
      EtGmPhd(const GmPhdCalibrations<4u, 2u> & calibrations)
        : BaseGmPhd<state_size, measurement_size, EtHypothesis<state_size>, EtInputHypothesis<state_size, measurement_size>>(calibrations) {}

      ~EtGmPhd(void) = default;

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

    protected:
      virtual std::pair<StateSizeVector, StateSizeMatrix> PredictHypothesis(const StateSizeVector & state, const StateSizeMatrix & covariance) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;
      
      void PredictBirths(void) {
        constexpr auto birth_objects_number = 100u;
        for (auto index = 0; index < birth_objects_number; index++) {
          Hypothesis birth_hypothesis;

          birth_hypothesis.weight = 2.0 / static_cast<double>(birth_objects_number);

          birth_hypothesis.state(0u) = pose_dist(e);
          birth_hypothesis.state(1u) = pose_dist(e);
          birth_hypothesis.state(2u) = velocity_dist(e);
          birth_hypothesis.state(3u) = velocity_dist(e);

          birth_hypothesis.covariance = 1.0 * StateSizeMatrix::Identity();

          const auto predicted_measurement = calibrations_.observation_matrix * birth_hypothesis.state;
          const auto innovation_covariance = calibrations_.measurement_covariance + calibrations_.observation_matrix * birth_hypothesis.covariance * calibrations_.observation_matrix.transpose();
          const auto kalman_gain = birth_hypothesis.covariance * calibrations_.observation_matrix.transpose() * innovation_covariance.inverse();
          const auto predicted_covariance = (StateSizeMatrix::Identity() - kalman_gain * calibrations_.observation_matrix) * birth_hypothesis.covariance;

          predicted_hypothesis_.push_back(birth_hypothesis);
        }
      }

      std::vector<Hypothesis> predicted_hypothesis_;

    private:
      void PredictExistingTargets(void) {
        // Prepare for prediction 
        PrepareTransitionMatrix();
        PrepareProcessNoiseMatrix();
        // Predict
        std::transform(hypothesis_.begin(), hypothesis_.end(),
          std::back_inserter(predicted_hypothesis_),
          [this](const Hypothesis & hypothesis) {
            // Predict kinematics
            const auto [predicted_state, predicted_state_covariance] = PredictHypothesis(hypothesis.state, hypothesis.covariance);
            const auto predicted_weight = calibrations_.ps * hypothesis.weight;
            const Hypothesis predicted_hypothesis = Hypothesis(predicted_weight, predicted_state, predicted_state_covariance);

            // Predict shape
            // ?

            return Hypothesis(predicted_state);
          }
        );
      }

      void Update(const std::vector<Measurement> & measurement) {
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

      std::vector<PredictedHypothesis> input_hypothesis_;
      DistancePartitioner<measurement_size> partitioner_;

      const double gamma_ = 1.0;
      const double lambda_ = 0.1;
      const double ck_ = 1.0e-9;
  };
} //  namespace mot

#endif  //  GM_PHD_INCLUDE_ET_GM_PHD_HPP_
