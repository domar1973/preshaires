#pragma once

#include <cmath>

namespace preshaires {

struct Energy {
  double eV;
};

struct MagneticField {
  double tesla;
};

struct Length {
  double meter;
};

struct Time {
  double second;
};

struct Rate {
  double per_second;
};

[[nodiscard]] constexpr Energy electron_volts(double value) {
  return Energy{value};
}

[[nodiscard]] constexpr Energy gigaelectron_volts(double value) {
  return Energy{value * 1.0e9};
}

[[nodiscard]] constexpr MagneticField tesla(double value) {
  return MagneticField{value};
}

[[nodiscard]] constexpr MagneticField nanotesla(double value) {
  return MagneticField{value * 1.0e-9};
}

[[nodiscard]] constexpr Length meters(double value) {
  return Length{value};
}

[[nodiscard]] constexpr Time seconds(double value) {
  return Time{value};
}

[[nodiscard]] inline bool is_finite(Energy value) {
  return std::isfinite(value.eV);
}

[[nodiscard]] inline bool is_finite(MagneticField value) {
  return std::isfinite(value.tesla);
}

[[nodiscard]] inline bool is_finite(Length value) {
  return std::isfinite(value.meter);
}

[[nodiscard]] inline bool is_finite(Time value) {
  return std::isfinite(value.second);
}

[[nodiscard]] inline bool is_finite(Rate value) {
  return std::isfinite(value.per_second);
}

}  // namespace preshaires
