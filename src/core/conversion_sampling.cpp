#include "preshaires/conversion_sampling.hpp"

#include <cmath>
#include <stdexcept>

namespace preshaires {
namespace {

void validate_localization_options(const LocalizationOptions& options) {
  if (!is_finite(options.absolute_path_tolerance) ||
      options.absolute_path_tolerance.meter < 0.0 ||
      !std::isfinite(options.relative_path_tolerance) ||
      options.relative_path_tolerance < 0.0) {
    throw std::invalid_argument(
        "Localization tolerances must be finite and nonnegative");
  }
  if (options.absolute_path_tolerance.meter == 0.0 &&
      options.relative_path_tolerance == 0.0) {
    throw std::invalid_argument(
        "At least one localization tolerance must be positive");
  }
  if (options.max_iterations == 0) {
    throw std::invalid_argument("Localization max_iterations must be positive");
  }
}

double checked_uniform_open01(RandomEngine& rng) {
  const double value = rng.uniform_open01();
  if (!std::isfinite(value) || value <= 0.0 || value >= 1.0) {
    throw std::runtime_error("RandomEngine returned a value outside (0, 1)");
  }
  return value;
}

bool length_tolerance_satisfied(double low, double high,
                                const LocalizationOptions& options) {
  const double width = high - low;
  const double tolerance =
      options.absolute_path_tolerance.meter +
      options.relative_path_tolerance * std::max(std::abs(low), std::abs(high));
  return width <= tolerance;
}

}  // namespace

ConversionSample sample_first_conversion(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, RandomEngine& rng,
    const IntegrationOptions& integration_options,
    const LocalizationOptions& localization_options) {
  validate_localization_options(localization_options);

  const double uniform = checked_uniform_open01(rng);
  const OpticalDepth target{-std::log(uniform)};
  const OpticalDepthResult total =
      optical_depth(photon_energy, trajectory, field, trajectory.length,
                    integration_options);
  const Probability total_probability =
      conversion_probability_from_optical_depth(total.optical_depth);
  std::size_t evaluations = total.evaluations;

  if (total.optical_depth.value == 0.0 ||
      target.value >= total.optical_depth.value) {
    return ConversionSample{false,
                            uniform,
                            target,
                            total.optical_depth,
                            total_probability,
                            std::nullopt,
                            std::nullopt,
                            evaluations};
  }

  double low = 0.0;
  double high = trajectory.length.meter;
  OpticalDepthResult high_tau = total;

  for (std::size_t iteration = 0;
       iteration < localization_options.max_iterations; ++iteration) {
    const double mid = 0.5 * (low + high);
    const OpticalDepthResult mid_tau =
        optical_depth(photon_energy, trajectory, field, meters(mid),
                      integration_options);
    evaluations += mid_tau.evaluations;

    if (mid_tau.optical_depth.value >= target.value) {
      high = mid;
      high_tau = mid_tau;
    } else {
      low = mid;
    }

    if (length_tolerance_satisfied(low, high, localization_options)) {
      const Length path = meters(high);
      return ConversionSample{true,
                              uniform,
                              target,
                              total.optical_depth,
                              total_probability,
                              path,
                              position_at(trajectory, path),
                              evaluations};
    }
  }

  throw std::runtime_error("Conversion localization exceeded max_iterations");
}

}  // namespace preshaires
