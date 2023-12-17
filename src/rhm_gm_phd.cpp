#include "rhm_gm_phd.hpp"

namespace mot {
  RhmGmPhd::RhmGmPhd() {};

  RhmGmPhd::~RhmGmPhd(void) = default;

  void RhmGmPhd::Run(const double timestamp, const std::vector<Measurement>& measurements) {
    time_delta_ = static_cast<float>(timestamp - prev_timestamp_);
    prev_timestamp_ = timestamp;

    MakePrediction();
    MakeCorrection(measurements);
    MakePruning();
    MakeMerging();
    ExtractObjects();
  }

  const RhmGmPhd::ObjectsList& RhmGmPhd::GetObjects(void) const {
    return objects_list_;
  }

  void RhmGmPhd::MakePrediction(void) {
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      std::back_inserter(predicted_hypothesis_list_),
      [this](const Hypothesis& hypothesis) -> Hypothesis {
        const auto predicted_weight = ps_ * hypothesis.weight;

        const KinematicStateSizeVector predicted_kinematic_state = transition_matrix_ * hypothesis.kinematic.value;
        const KinematicStateSizeMatrix predicted_kinematic_covariance = transition_matrix_ * hypothesis.kinematic.covariance * transition_matrix_.transpose() + q_;
        const ValueWithCovariance<kinematic_state_size> predicted_kinematic= {
          predicted_kinematic_state,      //  value
          predicted_kinematic_covariance  // kovariance
        };

        const ValueWithCovariance<extend_state_size> predicted_extend = {};

        Hypothesis predicted_hypothesis = {
          .weight = predicted_weight,     // weight
          .kinematic = predicted_kinematic,  // kinematic
          .extend = predicted_extend      // extend
        };
        
        return predicted_hypothesis;
      }
    );
  }

  void RhmGmPhd::MakeCorrection(const std::vector<Measurement>& measurements) {
    predicted_hypothesis_list_.clear();

    // Modify weights
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      std::back_inserter(predicted_hypothesis_list_),
      [this](const Hypothesis& hypothesis) -> Hypothesis {
        Hypothesis new_hypothesis = {
          .weight = 1.0f - (1.0f - std::exp(-gamma_)) * pd_ * hypothesis.weight,
          .kinematic = hypothesis.kinematic,
          .extend = hypothesis.extend
        };
        return new_hypothesis;
      }
    );

    // Go over partitioning
    auto l = 0u;
    for (const auto & partition : partitions_) {
      // for (const auto & point : partition.points) {
      for (auto j = 0; j < hypothesis_list_.size(); j++) {
        // TODO
        for (auto nz = 0u; nz < partition.points.size(); nz++) {

        }
      }
    };
  }

  void RhmGmPhd::MakePruning(void) {}

  void RhmGmPhd::MakeMerging(void) {}

  void RhmGmPhd::ExtractObjects(void) {
    objects_list_.clear();
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      std::back_inserter(objects_list_),
      [](const Hypothesis& hypothesis) -> Object {
        Object object;
        return object;
      }
    );
  };
}