#include "preshaires/conversion_probability.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "preshaires/constants.hpp"
#include "preshaires/pair_production.hpp"

namespace preshaires {
namespace {

struct Interval {
  double a;
  double b;
  double fa;
  double fm;
  double fb;
  double whole;
  double tolerance;
  int depth;
};

void validate_options(const IntegrationOptions& options) {
  if (!std::isfinite(options.relative_tolerance) ||
      !std::isfinite(options.absolute_tolerance) ||
      options.relative_tolerance < 0.0 || options.absolute_tolerance < 0.0) {
    throw std::invalid_argument("Integration tolerances must be finite and nonnegative");
  }
  if (options.relative_tolerance == 0.0 && options.absolute_tolerance == 0.0) {
    throw std::invalid_argument("At least one integration tolerance must be positive");
  }
  if (options.max_evaluations == 0) {
    throw std::invalid_argument("max_evaluations must be positive");
  }
}

void validate_trajectory(const StraightTrajectory& trajectory) {
  if (!is_finite(trajectory.length) || trajectory.length.meter < 0.0) {
    throw std::invalid_argument("Trajectory length must be finite and nonnegative");
  }
  (void)position_at(trajectory, meters(0.0));
}

double simpson(double a, double b, double fa, double fm, double fb) {
  return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

Probability probability_from_tau(double tau) {
  if (tau == 0.0) {
    return Probability{0.0};
  }
  const double probability = -std::expm1(-tau);
  return Probability{std::clamp(probability, 0.0, 1.0)};
}

}  // namespace

ConversionProbabilityResult conversion_probability(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, const IntegrationOptions& options) {
  if (!is_finite(photon_energy) || photon_energy.eV <= 0.0) {
    throw std::invalid_argument("Photon energy must be finite and positive");
  }
  validate_trajectory(trajectory);
  validate_options(options);

  const double length = trajectory.length.meter;
  std::size_t evaluations = 0;
  const auto evaluate = [&](double s) {
    if (evaluations >= options.max_evaluations) {
      throw std::runtime_error("Integration exceeded max_evaluations");
    }
    ++evaluations;
    const Length path = meters(s);
    const Position position = position_at(trajectory, path);
    const MagneticFieldVector local_field = field.field_at(position, path);
    const MagneticField b_perp =
        transverse_field(local_field, trajectory.direction);
    const Rate rate = erber1966_pair_production_rate(photon_energy, b_perp);
    if (!is_finite(rate) || rate.per_second < 0.0) {
      throw std::runtime_error("Nonfinite or negative pair-production rate");
    }
    const double integrand =
        rate.per_second / constants::speed_of_light_m_per_s;
    if (!std::isfinite(integrand) || integrand < 0.0) {
      throw std::runtime_error("Nonfinite or negative optical-depth integrand");
    }
    return integrand;
  };

  if (length == 0.0) {
    return ConversionProbabilityResult{OpticalDepth{0.0}, Probability{0.0}, 0};
  }

  const double fa = evaluate(0.0);
  const double fb = evaluate(length);
  const double midpoint = 0.5 * length;
  const double fm = evaluate(midpoint);
  const double whole = simpson(0.0, length, fa, fm, fb);
  const double requested_tolerance =
      options.absolute_tolerance + options.relative_tolerance * std::abs(whole);

  double integral = 0.0;
  std::vector<Interval> stack;
  stack.push_back(
      Interval{0.0, length, fa, fm, fb, whole, requested_tolerance, 0});

  constexpr int max_depth = 64;
  while (!stack.empty()) {
    const Interval current = stack.back();
    stack.pop_back();

    const double mid = 0.5 * (current.a + current.b);
    const double left_mid = 0.5 * (current.a + mid);
    const double right_mid = 0.5 * (mid + current.b);
    const double f_left_mid = evaluate(left_mid);
    const double f_right_mid = evaluate(right_mid);
    const double left =
        simpson(current.a, mid, current.fa, f_left_mid, current.fm);
    const double right =
        simpson(mid, current.b, current.fm, f_right_mid, current.fb);
    const double refined = left + right;
    const double error = std::abs(refined - current.whole) / 15.0;

    if (error <= current.tolerance) {
      integral += refined + (refined - current.whole) / 15.0;
      continue;
    }

    if (current.depth >= max_depth) {
      throw std::runtime_error("Adaptive Simpson integration exceeded max depth");
    }

    const double child_tolerance = 0.5 * current.tolerance;
    stack.push_back(Interval{mid, current.b, current.fm, f_right_mid,
                             current.fb, right, child_tolerance,
                             current.depth + 1});
    stack.push_back(Interval{current.a, mid, current.fa, f_left_mid,
                             current.fm, left, child_tolerance,
                             current.depth + 1});
  }

  if (!std::isfinite(integral) || integral < 0.0) {
    throw std::runtime_error("Integrated optical depth is invalid");
  }

  return ConversionProbabilityResult{OpticalDepth{integral},
                                     probability_from_tau(integral),
                                     evaluations};
}

}  // namespace preshaires
