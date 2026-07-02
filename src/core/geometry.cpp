#include "preshaires/geometry.hpp"

#include <cmath>
#include <stdexcept>

namespace preshaires {
namespace {

[[nodiscard]] bool finite3(double x, double y, double z) {
  return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

void validate_unit_direction(const Direction& direction) {
  if (!finite3(direction.x, direction.y, direction.z)) {
    throw std::invalid_argument("Direction components must be finite");
  }
  const double norm = std::hypot(direction.x, direction.y, direction.z);
  if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
    throw std::invalid_argument("Direction must be a unit vector");
  }
}

}  // namespace

Direction make_direction(double x, double y, double z) {
  if (!finite3(x, y, z)) {
    throw std::invalid_argument("Direction components must be finite");
  }

  const double norm = std::hypot(x, y, z);
  if (!std::isfinite(norm) || norm == 0.0) {
    throw std::invalid_argument("Direction norm must be finite and nonzero");
  }

  Direction direction{x / norm, y / norm, z / norm};
  const double normalized_norm =
      std::hypot(direction.x, direction.y, direction.z);
  if (!std::isfinite(normalized_norm) ||
      std::abs(normalized_norm - 1.0) > 1.0e-14) {
    throw std::runtime_error("Direction normalization failed");
  }
  return direction;
}

Position position_at(const StraightTrajectory& trajectory, Length path_length) {
  if (!is_finite(path_length) || path_length.meter < 0.0) {
    throw std::invalid_argument("Path length must be finite and nonnegative");
  }
  if (!is_finite(trajectory.length) || trajectory.length.meter < 0.0) {
    throw std::invalid_argument("Trajectory length must be finite and nonnegative");
  }
  if (!is_finite(trajectory.start.x) || !is_finite(trajectory.start.y) ||
      !is_finite(trajectory.start.z)) {
    throw std::invalid_argument("Trajectory start position must be finite");
  }
  validate_unit_direction(trajectory.direction);
  if (path_length.meter > trajectory.length.meter) {
    throw std::out_of_range("Path length is outside the trajectory");
  }

  return Position{
      meters(trajectory.start.x.meter + path_length.meter * trajectory.direction.x),
      meters(trajectory.start.y.meter + path_length.meter * trajectory.direction.y),
      meters(trajectory.start.z.meter + path_length.meter * trajectory.direction.z)};
}

MagneticField transverse_field(const MagneticFieldVector& field,
                               const Direction& direction) {
  if (!is_finite(field.x) || !is_finite(field.y) || !is_finite(field.z)) {
    throw std::invalid_argument("Magnetic field components must be finite");
  }
  validate_unit_direction(direction);

  const double b2 = field.x.tesla * field.x.tesla +
                    field.y.tesla * field.y.tesla +
                    field.z.tesla * field.z.tesla;
  const double b_dot_n = field.x.tesla * direction.x +
                         field.y.tesla * direction.y +
                         field.z.tesla * direction.z;
  double b_perp2 = b2 - b_dot_n * b_dot_n;
  if (b_perp2 < 0.0 && b_perp2 > -1.0e-24 * (1.0 + b2)) {
    b_perp2 = 0.0;
  }
  if (b_perp2 < 0.0) {
    throw std::runtime_error("Computed negative transverse field norm");
  }

  return tesla(std::sqrt(b_perp2));
}

}  // namespace preshaires
