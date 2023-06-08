#ifndef GM_PHD_INCLUDE_BASE_GM_PHD_HPP_
#define GM_PHD_INCLUDE_BASE_GM_PHD_HPP_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <random>
#include <ranges>
#include <vector>

#include <Eigen/Dense>

#include "gm_phd_calibrations.hpp"
#include "value_with_covariance.hpp"

namespace mot {
  std::random_device r;
  std::default_random_engine e(r());

  std::uniform_real_distribution<double> pose_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> velocity_dist(-1.0, 1.0);

  template <size_t state_size, size_t measurement_size, typename Hypothesis, typename PredictedHypothesis>
  class BaseGmPhd {
    public:
      using StateSizeVector = Eigen::Vector<double, state_size>;
      using StateSizeMatrix = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementSizeVector = Eigen::Vector<double, measurement_size>;
      using MeasurementSizeMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

      using Object = ValueWithCovariance<state_size>;
      using Measurement = ValueWithCovariance<measurement_size>;

    public:
      explicit BaseGmPhd(const GmPhdCalibrations<state_size, measurement_size> & calibrations)
        : calibrations_{calibrations} {}

      virtual ~BaseGmPhd(void) = default;

      virtual void Run(const double timestamp, const std::vector<Measurement> & measurements) = 0;

      const std::vector<Object> & GetObjects(void) const {
        return objects_;
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
      virtual std::pair<StateSizeVector, StateSizeMatrix> PredictHypothesis(const StateSizeVector & state, const StateSizeMatrix & covariance) = 0;
      virtual void PrepareTransitionMatrix(void) = 0;
      virtual void PrepareProcessNoiseMatrix(void) = 0;
      virtual void PredictBirths(void) = 0;

      std::vector<PredictedHypothesis> predicted_hypothesis_;

      double time_delta = 0.0;
      GmPhdCalibrations<state_size, measurement_size> calibrations_;
      StateSizeMatrix transition_matrix_ = StateSizeMatrix::Zero();
      StateSizeMatrix process_noise_covariance_matrix_ = StateSizeMatrix::Zero();

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

      virtual void PredictExistingTargets(void) = 0;

      virtual void Update(const std::vector<Measurement> & measurements) = 0;

      virtual void UpdateExistedHypothesis(void) = 0;

      virtual void MakeMeasurementUpdate(const std::vector<Measurement> & measurements) = 0;

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

          merged_hypothesis.push_back(Hypothesis(merged_weight, merged_state, merged_covariance));
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
            object.value = hypothesis.state;
            object.covariance = hypothesis.covariance;
            objects_.push_back(object);
          }
        }
      }

      double NormPdf(const MeasurementSizeVector & z, const MeasurementSizeVector & nu, const MeasurementSizeMatrix & cov) const {
        const auto diff = z - nu;
        const auto c = 1.0 / (std::sqrt(std::pow(std::numbers::pi, measurement_size) * cov.determinant()));
        const auto e = std::exp(-0.5 * diff.transpose() * cov.inverse() * diff);
        return c * e;
      }

      double prev_timestamp_ = 0.0;
      bool is_initialized_ = false;
      std::vector<Object> objects_;
      std::vector<Hypothesis> hypothesis_;
  };
};  //  namespace mot

#endif  //  GM_PHD_INCLUDE_BASE_GM_PHD_HPP_
