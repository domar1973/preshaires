#pragma once

#include <string_view>

#include "preshaires/units.hpp"

namespace preshaires {

enum class PairProductionModel {
  Erber1966,
};

struct PhysicsMetadata {
  PairProductionModel model;
  std::string_view reference;
  std::string_view equations;
  std::string_view assumptions;
  std::string_view validity_range;
};

[[nodiscard]] PhysicsMetadata erber1966_metadata();

[[nodiscard]] double photon_chi(Energy photon_energy,
                                MagneticField transverse_magnetic_field);

[[nodiscard]] Rate erber1966_pair_production_rate(
    Energy photon_energy, MagneticField transverse_magnetic_field);

}  // namespace preshaires
