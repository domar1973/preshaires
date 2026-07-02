#pragma once

#include <cstddef>
#include <optional>

#include "preshaires/conversion_probability.hpp"
#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
#include "preshaires/random.hpp"
#include "preshaires/units.hpp"

namespace preshaires {

struct LocalizationOptions {
  Length absolute_path_tolerance;
  double relative_path_tolerance;
  std::size_t max_iterations;
};

struct ConversionSample {
  bool converted;
  double random_uniform;
  OpticalDepth target_optical_depth;
  OpticalDepth total_optical_depth;
  Probability total_probability;
  std::optional<Length> path_length;
  std::optional<Position> position;
  std::size_t evaluations;
};

[[nodiscard]] ConversionSample sample_first_conversion(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, RandomEngine& rng,
    const IntegrationOptions& integration_options,
    const LocalizationOptions& localization_options);

}  // namespace preshaires
