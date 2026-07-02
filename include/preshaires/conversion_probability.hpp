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

struct OpticalDepthResult {
  OpticalDepth optical_depth;
  std::size_t evaluations;
};

[[nodiscard]] Probability conversion_probability_from_optical_depth(
    OpticalDepth optical_depth);

[[nodiscard]] OpticalDepthResult optical_depth(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, Length upper_path_length,
    const IntegrationOptions& options);

[[nodiscard]] ConversionProbabilityResult conversion_probability(
    Energy photon_energy, const StraightTrajectory& trajectory,
    const MagneticFieldModel& field, const IntegrationOptions& options);

}  // namespace preshaires
