#include "rhm_gm_phd.hpp"

#include <random>

namespace mot {
  RhmGmPhd::RhmGmPhd() {};

  RhmGmPhd::~RhmGmPhd(void) = default;

  void RhmGmPhd::Run(const double timestamp, const std::vector<Measurement>& measurements) {
    time_delta_ = static_cast<float>(timestamp - prev_timestamp_);
    prev_timestamp_ = timestamp;

    predicted_hypothesis_list_.clear();

    MakeBirth();
    MakePrediction();
    MakeCorrection(measurements);
    MakePruning();
    MakeMerging();
    ExtractObjects();
  }

  const RhmGmPhd::ObjectsList& RhmGmPhd::GetObjects(void) const {
    return objects_list_;
  }

  void RhmGmPhd::MakeBirth(void) {
    static std::random_device r;
    static std::default_random_engine e(r());

    static std::uniform_real_distribution<double> pose_dist(-100.0, 100.0);
    static std::uniform_real_distribution<double> velocity_dist(-5.0, 5.0);

    constexpr auto birth_objects_number = 100u;
    for (auto index = 0; index < birth_objects_number; index++) {
      Hypothesis birth_hypothesis;

      birth_hypothesis.weight = 2.0 / static_cast<double>(birth_objects_number);

      birth_hypothesis.kinematic.value(0u) = pose_dist(e);
      birth_hypothesis.kinematic.value(1u) = pose_dist(e);
      birth_hypothesis.kinematic.value(2u) = velocity_dist(e);
      birth_hypothesis.kinematic.value(3u) = velocity_dist(e);

      birth_hypothesis.kinematic.covariance = KinematicStateSizeMatrix::Identity();

      birth_hypothesis.extend.value(0u) = 0.5f;
      birth_hypothesis.extend.value(1u) = 0.5f;
      birth_hypothesis.extend.value(2u) = 0.5f;

      birth_hypothesis.extend.covariance = ExtendStateSizeMatrix::Identity();

      predicted_hypothesis_list_.push_back(birth_hypothesis);
    }
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

  void RhmGmPhd::MakePruning(void) {
    std::vector<Hypothesis> pruned_hypothesis;
    std::copy_if(predicted_hypothesis_list_.begin(), predicted_hypothesis_list_.end(),
      std::back_inserter(pruned_hypothesis),
      [this](const Hypothesis & hypothesis) -> bool {
        return hypothesis.weight >= 0.25; // TODO: add calibration calibrations_.truncation_threshold;
      }
    );

    predicted_hypothesis_list_.resize(pruned_hypothesis.size());
    std::copy(pruned_hypothesis.begin(), pruned_hypothesis.end(), predicted_hypothesis_list_.begin());
  }

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