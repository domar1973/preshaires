#pragma once

#include <cstddef>

#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
#include "preshaires/units.hpp"

namespace preshaires {

struct OpticalDepth {
  double value;
};

struct Probability {
  double value;
};

struct IntegrationOptions {
  double relative_tolerance;
  double absolute_tolerance;
  std::size_t max_evaluations;
};

struct ConversionProbabilityResult {
  OpticalDepth optical_depth;
  Probability probability;
  std::size_t evaluations;
};

[[nodiscard]] ConversionProbabilityResult conversion_probability(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, const IntegrationOptions& options);

}  // namespace preshaires
