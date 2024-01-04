#include "rhm_gm_phd.hpp"

#include <iterator>
#include <random>

namespace mot {
  constexpr auto pi = 3.141592653589793;

  RhmGmPhd::RhmGmPhd(const RhmGmPhdCalibrations& calibrations)
    : calibrations_{calibrations} {};

  RhmGmPhd::~RhmGmPhd(void) = default;

  void RhmGmPhd::Run(const double timestamp, const std::vector<Measurement>& measurements) {
    time_delta_ = static_cast<double>(timestamp - prev_timestamp_);
    prev_timestamp_ = timestamp;

    predicted_hypothesis_list_.clear();

    MakeBirth();
    MakePrediction();
    MakePartitioning(measurements);
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

    static std::uniform_real_distribution<double> pose_dist(-300.0, 300.0);
    static std::uniform_real_distribution<double> velocity_dist(-5.0, 5.0);

    constexpr auto birth_objects_number = 1000u;
    for (auto index = 0; index < birth_objects_number; index++) {
      Hypothesis birth_hypothesis;

      birth_hypothesis.weight = 0.2;//2.0 / static_cast<double>(birth_objects_number);

      birth_hypothesis.kinematic.value(0u) = pose_dist(e);
      birth_hypothesis.kinematic.value(1u) = pose_dist(e);
      birth_hypothesis.kinematic.value(2u) = velocity_dist(e);
      birth_hypothesis.kinematic.value(3u) = velocity_dist(e);

      birth_hypothesis.kinematic.covariance = KinematicStateSizeMatrix::Identity();

      birth_hypothesis.extend.value(0u) = 50.0;
      birth_hypothesis.extend.value(1u) = 50.0;
      birth_hypothesis.extend.value(2u) = 0.0;

      birth_hypothesis.extend.covariance = ExtendStateSizeMatrix::Identity();

      birth_hypothesis.s = 0.7;
      birth_hypothesis.var_s = 0.08;

      birth_hypothesis.var_zero = 1.0;

      hypothesis_list_.push_back(birth_hypothesis);
    }
  }

  void RhmGmPhd::MakePrediction(void) {
    predicted_hypothesis_list_.clear();
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      std::back_inserter(predicted_hypothesis_list_),
      [this](const Hypothesis& hypothesis) -> Hypothesis {
        const auto predicted_weight = calibrations_.ps * hypothesis.weight;

        const KinematicStateSizeVector predicted_kinematic_state = transition_matrix_ * hypothesis.kinematic.value;
        const KinematicStateSizeMatrix predicted_kinematic_covariance = transition_matrix_ * hypothesis.kinematic.covariance * transition_matrix_.transpose() + q_;
        const ValueWithCovariance<kinematic_state_size> predicted_kinematic = {
          predicted_kinematic_state,      //  value
          predicted_kinematic_covariance  // covariance
        };

        const ValueWithCovariance<extend_state_size> predicted_extend = hypothesis.extend;

        Hypothesis predicted_hypothesis = {
          .weight = predicted_weight,       // weight
          .kinematic = predicted_kinematic, // kinematic
          .extend = predicted_extend,       // extend
          .s = hypothesis.s,
          .var_s = hypothesis.var_s,
          .var_zero = hypothesis.var_zero
        };
        
        return predicted_hypothesis;
      }
    );
    jk_ = predicted_hypothesis_list_.size();
  }

  void RhmGmPhd::MakePartitioning(const std::vector<Measurement>& measurements) {
    partitions_ = distance_partitioner_.MakePartitioning(measurements);
  }

  void RhmGmPhd::MakeCorrection(const std::vector<Measurement>& measurements) {
    // Modify weights
    std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
      predicted_hypothesis_list_.begin(),
      [this](const Hypothesis& hypothesis) -> Hypothesis {
        Hypothesis new_hypothesis = {
          .weight = 1.0 - (1.0 - std::exp(-gamma_)) * calibrations_.pd * hypothesis.weight,
          .kinematic = hypothesis.kinematic,
          .extend = hypothesis.extend,
          .s = hypothesis.s,
          .var_s = hypothesis.var_s,
          .var_zero = hypothesis.var_zero
        };
        return new_hypothesis;
      }
    );

    // Go over partitioning
    std::copy(predicted_hypothesis_list_.begin(), predicted_hypothesis_list_.end(), std::back_inserter(hypothesis_list_));
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
        scale_ = hypothesis_list_.at(j).s;
        scale_var_ = hypothesis_list_.at(j).var_s;

        for (auto nz = 0u; nz < partition.points.size(); nz++) {
          const auto measurement = measurements.at(partition.points.at(nz));
          MakeKinematicCorrection(measurement);
          MakeExtendCorrection(hypothesis_list_.at(j), measurement);
        }

        // Line 24, 25
        pl_ *= std::exp(-gamma_) * calibrations_.pd * gamma_;
        Hypothesis new_hypothesis;
        new_hypothesis.weight = pl_;
        new_hypothesis.kinematic.value = mu_.tail<kinematic_state_size>();
        new_hypothesis.kinematic.covariance = px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size);
        new_hypothesis.extend.value = mu_.head<extend_state_size>();
        new_hypothesis.extend.covariance = px_.block<extend_state_size, extend_state_size>(0u, 0u);
        new_hypothesis.s = 0.7;
        new_hypothesis.var_s = 0.08;

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

  void RhmGmPhd::MakeKinematicCorrection(const Measurement& measurement) {
    const Eigen::Vector<double, 4u> kinematic_state = mu_.tail<kinematic_state_size>();
    const Eigen::Matrix<double, 4u, 4u>  kinematic_state_covariance = px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size);

    // Make correction
    const auto& h = calibrations_.observation_matrix;
    const Eigen::Vector<double, 2u> innovation = measurement.value - h * kinematic_state;
    const Eigen::Matrix<double, 2u, 2u> innovation_covariance = h * kinematic_state_covariance * h.transpose() + measurement.covariance;
    const Eigen::Matrix<double, 4u, 2u> kalman_gain = kinematic_state_covariance * h.transpose() * innovation_covariance.inverse();

    mu_.tail<kinematic_state_size>() += kalman_gain * innovation;
    px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size) = (KinematicStateSizeMatrix::Identity() - kalman_gain * h) * kinematic_state_covariance;

    // Force symetricity
    //px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size) = MakeMatrixSymetric<kinematic_state_size>(px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size));
  }

  void RhmGmPhd::MakeExtendCorrection(const Hypothesis& hypothesis, const Measurement& measurement) {
    // Line 12
    // Eigen::Vector<double, augmented_size> mu_a = Eigen::Vector<double, augmented_size>::Zero();
    // mu_a.head<nx>() = mu_;
    // mu_a(nx) = hypothesis.s;//s?
    // Eigen::Matrix<double, augmented_size, augmented_size> c_a = Eigen::Matrix<double, augmented_size, augmented_size>::Identity();
    // c_a.block<nx, nx>(0u, 0u) = px_;
    // c_a(nx, nx) = hypothesis.var_s;// sigma_s
    // c_a.block<measurement_size, measurement_size>(nx + 1u, nx + 1u) = measurement.covariance;
    // // Line 13
    // const auto theta = std::atan2(measurement.value(1) - m_(1), measurement.value(0) - m_(0));
    // // Line 14, 15
    // const auto ut_points = GetUkfPointsAndWeights(mu_a, c_a, theta, measurement);
    // // Line 16, 17, 18, 19
    // const auto k = UtCorrectionStep(ut_points);
    // // line 20
    // mu_a += k * (0.0f - z_mean_);
    // c_a -= k * s_ * k.transpose();

    // px_.block<extend_state_size, extend_state_size>(0u, 0u) = c_a.block<extend_state_size, extend_state_size>(0u, 0u);
    // px_.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size) = c_a.block<kinematic_state_size, kinematic_state_size>(extend_state_size, extend_state_size);
    // mu_.head<extend_state_size>() = mu_a.head<extend_state_size>();
    // mu_.tail<kinematic_state_size>() = mu_a.tail<kinematic_state_size>();

    // pl_ *= std::max((1.0 / (std::sqrt(s_ * 2.0 * pi))) * std::exp(-0.5 * z_mean_ * s_ * z_mean_), 0.01);




    /* Noise */
    static Eigen::Vector3d measurement_noise_mean = Eigen::Vector3d::Zero();
    measurement_noise_mean(0u) = scale_;

    static Eigen::Matrix3d measurement_noise_cov = Eigen::Matrix3d::Zero();
    measurement_noise_cov(0u, 0u) = scale_var_;
    measurement_noise_cov.block<2u, 2u>(1u, 1u) = measurement.covariance;

    /* Prepare current extended state */
    Eigen::Vector<double, extend_state_size + 3u> p_ukf;
    p_ukf.head(extend_state_size) = mu_.head<extend_state_size>();
    p_ukf.tail(3u) = measurement_noise_mean;

    Eigen::Matrix<double, extend_state_size + 3u, extend_state_size + 3u> c_ukf;
    c_ukf.block<extend_state_size, extend_state_size>(0u, 0u) = px_.block<extend_state_size, extend_state_size>(0u, 0u);
    c_ukf.block<3u, 3u>(extend_state_size, extend_state_size) = measurement_noise_cov;

    // Line 13
    const auto theta = std::atan2(measurement.value(1) - m_(1), measurement.value(0) - m_(0));
    // Line 14, 15
    const auto ut_points = GetUkfPointsAndWeights(p_ukf, c_ukf, theta, measurement);
    // Line 16, 17, 18, 19
    const auto k = UtCorrectionStep(ut_points);

    // line 20
    mu_.head<extend_state_size>() += k.head<extend_state_size>() * (0.0f - z_mean_);
    px_.block<extend_state_size, extend_state_size>(0u, 0u) -= k.head<extend_state_size>() * z_mean_var_ * k.head<extend_state_size>().transpose();

    // scale_ += k.tail<1u>().norm() * (0.0f - z_mean_);
    // scale_var_ -= k.tail<1u>().norm() * z_mean_var_ * k.tail<1u>().transpose().norm();

    pl_ *= std::max((1.0 / (std::sqrt(z_mean_var_ * 2.0 * pi))) * std::exp(-0.5 * z_mean_ * z_mean_var_ * z_mean_), 0.01);
  }

  RhmGmPhd::UtPoints RhmGmPhd::GetUkfPointsAndWeights(const UkfVector& mu_a, const UkfMatrix& c_a, const double theta, const Measurement& measurement) {
    UtPoints ut_points;
    const auto n = extend_state_size + 3u;
    const auto lambda = std::pow(alpha_, 2u) + (static_cast<double>(n) + kappa_)  -  static_cast<double>(n);

    ut_points.points_number = 2u * n + 1u;

    // Line 14
    ut_points.wm.resize(2u * n + 1u);
    ut_points.wm.at(0u) = lambda / (static_cast<double>(n) + lambda);
    std::fill_n(std::next(ut_points.wm.begin(), 1), 2u * n, 1.0 / (2.0 * (static_cast<double>(n)  + lambda)));

    ut_points.wc.resize(2u * n + 1u);
    ut_points.wc.at(0u) = (lambda / (static_cast<double>(n) + lambda)) + (1 - std::pow(alpha_, 2.0) + beta_);
    std::fill_n(std::next(ut_points.wc.begin(), 1), 2u * n, 1.0 / (2.0 * (static_cast<double>(n)  + lambda)));

    const UkfMatrix chol = c_a.llt().matrixL();
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
      [theta, measurement, this](const UkfVector& mu_a) -> double {
        const auto a = mu_a(0u);
        const auto b = mu_a(1u);
        const auto phi = mu_a(2u);
        const auto s = mu_a(3u);
        const auto v = mu_a.tail<2u>();
        Eigen::Vector2d e;
        e << std::cos(theta), std::sin(theta);
        Eigen::Vector2d m;
        m << mu_(extend_state_size), mu_(extend_state_size + 1u);

        const auto r = a * b / std::sqrt(std::pow(a * std::sin(theta - phi), 2.0) + std::pow(b * std::cos(theta - phi), 2.0));

        double expected_pseudomeasurement = std::pow(s * r, 2)
          + 2.0 * s * r * (e.transpose() * v).norm()
          + std::pow(v.norm(), 2)
          - std::pow((measurement.value - m).norm(), 2);

        return expected_pseudomeasurement;
      }
    );

    return ut_points;
  }

  RhmGmPhd::UkfVector RhmGmPhd::UtCorrectionStep(const RhmGmPhd::UtPoints& ut_points) {
    z_mean_ = 0.0;
    for (auto index = 0u; index < ut_points.points_number; index++)
      z_mean_ += ut_points.wm.at(index) * ut_points.z_pts.at(index);

    z_mean_var_ = 0.0;
    UkfVector p = UkfVector::Zero();
    for (auto index = 0u; index < ut_points.points_number; index++) {
      z_mean_var_ += ut_points.wc.at(index) * std::pow(ut_points.z_pts.at(index), 2);
      p += ut_points.wc.at(index) * ut_points.x_pts.at(index) * ut_points.z_pts.at(index);
    }

    const UkfVector k = p / z_mean_var_;
    return k;
  }

  void RhmGmPhd::MakePruning(void) {
    pruned_hypothesis_.clear();
    std::copy_if(predicted_hypothesis_list_.begin(), predicted_hypothesis_list_.end(),
      std::back_inserter(pruned_hypothesis_),
      [this](const Hypothesis& hypothesis) -> bool {
        return hypothesis.weight > calibrations_.truncation_threshold;
      }
    );
  }

  void RhmGmPhd::MakeMerging(void) {
    hypothesis_list_.resize(pruned_hypothesis_.size());
    std::copy(pruned_hypothesis_.begin(), pruned_hypothesis_.end(), hypothesis_list_.begin());
  }

  void RhmGmPhd::ExtractObjects(void) {
    objects_list_.clear();
    //objects_list_.resize(hypothesis_list_.size());
    for (const auto hypothesis : hypothesis_list_) {
      Object object = {
        .kinematic = hypothesis.kinematic,
        .extend = hypothesis.extend
      };
      objects_list_.push_back(object);
    }

    // std::transform(hypothesis_list_.begin(), hypothesis_list_.end(),
    //   std::back_inserter(objects_list_),
    //   [](const Hypothesis& hypothesis) -> Object {
    //     Object object = {
    //       .kinematic = hypothesis.kinematic,
    //       .extend = hypothesis.extend
    //     };
    //     return object;
    //   }
    // );
  };
} //  namespace mot
