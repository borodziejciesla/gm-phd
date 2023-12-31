#include "rhm_gm_phd.hpp"

#include <iterator>
#include <random>

namespace mot {
  constexpr auto pi = 3.141592653589793;

  RhmGmPhd::RhmGmPhd() {};

  RhmGmPhd::~RhmGmPhd(void) = default;

  void RhmGmPhd::Run(const double timestamp, const std::vector<Measurement>& measurements) {
    time_delta_ = static_cast<double>(timestamp - prev_timestamp_);
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

      birth_hypothesis.extend.value(0u) = 0.5;
      birth_hypothesis.extend.value(1u) = 0.5;
      birth_hypothesis.extend.value(2u) = 0.5;

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
        const ValueWithCovariance<kinematic_state_size> predicted_kinematic = {
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
    jk_ = predicted_hypothesis_list_.size();
  }

  void RhmGmPhd::MakePartitioning(void) {
    //
  }

  void RhmGmPhd::MakeCorrection(const std::vector<Measurement>& measurements) {
    predicted_hypothesis_list_.clear();

    // Modify weights
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      std::back_inserter(predicted_hypothesis_list_),
      [this](const Hypothesis& hypothesis) -> Hypothesis {
        Hypothesis new_hypothesis = {
          .weight = 1.0 - (1.0 - std::exp(-gamma_)) * pd_ * hypothesis.weight,
          .kinematic = hypothesis.kinematic,
          .extend = hypothesis.extend
        };
        return new_hypothesis;
      }
    );

    // Go over partitioning
    auto l = 0u;
    for (const auto & partition : partitions_) {
      l++;
      for (auto j = 0; j < hypothesis_list_.size(); j++) {
        // Line 10
        px_.block<extend_state_size, extend_state_size>(0u, 0u) = hypothesis_list_.at(j).extend.covariance;
        px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size) = hypothesis_list_.at(j).kinematic.covariance;
        mu_.head<extend_state_size>() = hypothesis_list_.at(j).extend.value;
        mu_.tail<kinematic_state_size>() = hypothesis_list_.at(j).kinematic.value;
        m_ << mu_(3u), mu_(4u);
        pl_ = 1.0;

        for (auto nz = 0u; nz < partition.points.size(); nz++) {
          ProcessMeasurementWithHypothesis(hypothesis_list_.at(j), measurements.at(partition.points.at(nz)));
        }

        // Line 24, 25
        pl_ *= std::exp(-gamma_) * pd_ * gamma_;
        Hypothesis new_hypothesis;
        new_hypothesis.weight = pl_;
        new_hypothesis.kinematic.value = mu_.tail<kinematic_state_size>();
        new_hypothesis.kinematic.covariance = px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size);
        new_hypothesis.extend.value = mu_.head<extend_state_size>();
        new_hypothesis.extend.covariance = px_.block<extend_state_size, extend_state_size>(0u, 0u);

        predicted_hypothesis_list_.push_back(new_hypothesis);
      }
      // Line 27
      auto dw = (partition.points.size() == 1) ? 1.0 : 0.0;
      for (auto j = 0u; j < jk_; j++) {
        dw += predicted_hypothesis_list_.at(j + jk_ * l).weight;
      }

      for (auto j = 0u; j < jk_; j++) {
        predicted_hypothesis_list_.at(j + jk_ * l).weight /= dw;
      }
    };
  }

  void RhmGmPhd::ProcessMeasurementWithHypothesis(const Hypothesis& hypothesis, const Measurement& measurement) {
    // Line 12
    Eigen::Vector<double, augmented_size> mu_a = Eigen::Vector<double, augmented_size>::Zero();
    mu_a.head<nx>() = mu_;
    mu_a(nx) = 0.0;//s?
    Eigen::Matrix<double, augmented_size, augmented_size> c_a = Eigen::Matrix<double, augmented_size, augmented_size>::Identity();
    c_a.block<nx, nx>(0u, 0u) = px_;
    c_a(nx, nx) = 1.0;// sigma_s
    c_a.block<measurement_size, measurement_size>(nx + 1u, nx + 1u) = measurement.covariance;
    // Line 13
    const Eigen::Vector2d diff = measurement.value - m_;
    const auto dx = diff(0u) * std::cos(mu_(3u)) - diff(1u) * std::sin(mu_(3u));
    const auto dy = diff(0u) * std::sin(mu_(3u)) + diff(1u) * std::cos(mu_(3u));
    const auto theta = std::atan2(dy, dx);
    // Line 14, 15
    const auto ut_points = GetUkfPointsAndWeights(mu_a, c_a, theta, measurement);
    // Line 16, 17, 18, 19
    const auto k = UtCorrectionStep(ut_points);
    // line 20
    mu_a += k * (0.0f - z_mean_);
    c_a -= k * s_ * k.transpose();

    px_.block<extend_state_size, extend_state_size>(0u, 0u) = c_a.block<extend_state_size, extend_state_size>(0u, 0u);
    px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size) = c_a.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size);
    mu_.head<extend_state_size>() = mu_a.head<extend_state_size>();
    mu_.tail<kinematic_state_size>() = mu_a.tail<kinematic_state_size>();

    pl_ *= (1.0 / (std::sqrt(s_ * 2.0 * pi))) * std::exp(-0.5 * z_mean_ * s_ * z_mean_);
  }

  RhmGmPhd::UtPoints RhmGmPhd::GetUkfPointsAndWeights(const AugmentedVector& mu_a, const AugmentedMatrix& c_a, const double theta, const Measurement& measurement) {
    UtPoints ut_points;
    const auto n = 5u;
    const auto lambda = std::pow(alpha_, 2u) + (static_cast<double>(n) + kappa_)  -  static_cast<double>(n);

    // Line 14
    ut_points.wm.resize(2u * n + 1u);
    ut_points.wm.at(0u) = lambda / (static_cast<double>(n) + lambda);
    std::fill_n(std::next(ut_points.wm.begin(), 1), 2u * n, 1.0 / (2.0 * (static_cast<double>(n)  + lambda)));

    ut_points.wc.resize(2u * n + 1u);
    ut_points.wc.at(0u) = (lambda / (static_cast<double>(n) + lambda)) + (1 - std::pow(alpha_, 2.0) + beta_);
    std::fill_n(std::next(ut_points.wc.begin(), 1), 2u * n, 1.0 / (2.0 * (static_cast<double>(n)  + lambda)));

    const AugmentedMatrix chol = c_a.llt().matrixL();
    const auto sigma_points_cholesky_part = std::sqrt(static_cast<double>(n) + lambda) * chol;

    ut_points.x_pts.resize(2u * n + 1u);
    ut_points.x_pts.at(0) = mu_a;
    for (auto index = 0u; index < n; index++) {
      const auto col_index = sigma_points_cholesky_part.col(index);
      ut_points.x_pts.at(index + 1u) = mu_a + col_index;
      ut_points.x_pts.at(n + index + 1u) = mu_a - col_index;
    }

    // Line 15
    std::transform(ut_points.x_pts.begin(), ut_points.x_pts.end(),
      std::back_inserter(ut_points.z_pts),
      [theta, measurement, this](const Eigen::Vector<double, augmented_size>& mu_a) -> double {
        const auto a = mu_a(0u);
        const auto b = mu_a(1u);
        const auto phi = mu_a(2u);
        const auto s = mu_a(augmented_size);
        const auto v = mu_a.tail<2u>();
        Eigen::Vector2d e;
        e << std::cos(theta), std::sin(theta);
        Eigen::Vector2d m;
        m << mu_a(kinematic_state_size), mu_a(kinematic_state_size + 1u);

        const auto r = a * b / std::sqrt(std::pow(a * std::sin(theta - phi), 2.0) + std::pow(b * std::cos(theta - phi), 2.0));

        double expected_pseudomeasurement = std::pow(s * r, 2)
          + 2.0 * s * r * (e.transpose() * v).norm()
          + v.norm()
          + (measurement.value - m).norm();

        return expected_pseudomeasurement;
      }
    );
  }

  RhmGmPhd::AugmentedVector RhmGmPhd::UtCorrectionStep(const RhmGmPhd::UtPoints& ut_points) {
    z_mean_ = 0.0;
    for (auto index = 0u; index < ut_points.points_number; index++)
      z_mean_ += ut_points.wm.at(index) * ut_points.z_pts.at(index);

    s_ = 0.0;
    AugmentedVector p = AugmentedVector::Zero();
    for (auto index = 0u; index < ut_points.points_number; index++) {
      s_ += ut_points.wc.at(index) * std::pow(ut_points.z_pts.at(index), 2);
      p += ut_points.wc.at(index) * ut_points.x_pts.at(index) * ut_points.z_pts.at(index);
    }

    const auto k = p / s_;
    return k;
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