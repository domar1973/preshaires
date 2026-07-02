#pragma once

#include "preshaires/units.hpp"

namespace preshaires {

struct Position {
  Length x;
  Length y;
  Length z;
};

struct Direction {
  double x;
  double y;
  double z;
};

struct MagneticFieldVector {
  MagneticField x;
  MagneticField y;
  MagneticField z;
};

struct StraightTrajectory {
  Position start;
  Direction direction;
  Length length;
};

[[nodiscard]] Direction make_direction(double x, double y, double z);

[[nodiscard]] Position position_at(const StraightTrajectory& trajectory,
                                   Length path_length);

[[nodiscard]] MagneticField transverse_field(
    const MagneticFieldVector& field, const Direction& direction);

}  // namespace preshaires
