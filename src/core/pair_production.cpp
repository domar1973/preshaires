#include "preshaires/pair_production.hpp"

#include <cmath>
#include <limits>

#include "preshaires/constants.hpp"

namespace preshaires {
namespace {

[[nodiscard]] double electron_rest_energy_joule() {
  return constants::electron_rest_energy.eV * constants::joule_per_electron_volt;
}

[[nodiscard]] bool invalid_nonpositive(Energy energy) {
  return !is_finite(energy) || energy.eV <= 0.0;
}

[[nodiscard]] bool invalid_field(MagneticField field) {
  return !is_finite(field);
}

}  // namespace

PhysicsMetadata erber1966_metadata() {
  return PhysicsMetadata{
      PairProductionModel::Erber1966,
      "T. Erber, Rev. Mod. Phys. 38, 626 (1966), Sec. 3B, Eq. (3.4).",
      "chi = 0.5 * (E_gamma / (m_e c^2)) * (B_perp / B_crit); "
      "Gamma = 0.16 * alpha * m_e c^2 / hbar * "
      "(m_e c^2 / E_gamma) * K_{1/3}(2/(3 chi))^2.",
      "Uniform local transverse magnetic field, H << H_crit, ultrarelativistic "
      "photon, leading Erber approximation, no bound-state resonances.",
      "Intended for chi in the asymptotic Erber approximation regime and "
      "geomagnetic fields far below B_crit; validate before extrapolating."};
}

double photon_chi(Energy photon_energy,
                  MagneticField transverse_magnetic_field) {
  if (invalid_nonpositive(photon_energy) || invalid_field(transverse_magnetic_field)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (transverse_magnetic_field.tesla == 0.0) {
    return 0.0;
  }

  const double energy_ratio =
      photon_energy.eV / constants::electron_rest_energy.eV;
  const double field_ratio =
      std::abs(transverse_magnetic_field.tesla) /
      constants::critical_magnetic_field.tesla;
  return 0.5 * energy_ratio * field_ratio;
}

Rate erber1966_pair_production_rate(
    Energy photon_energy, MagneticField transverse_magnetic_field) {
  const double chi = photon_chi(photon_energy, transverse_magnetic_field);
  if (!std::isfinite(chi)) {
    return Rate{std::numeric_limits<double>::quiet_NaN()};
  }
  if (chi <= 0.0) {
    return Rate{0.0};
  }

  constexpr double coefficient = 0.16;
  const double prefactor =
      coefficient * constants::fine_structure_constant *
      electron_rest_energy_joule() /
      constants::reduced_planck_constant_J_s;
  const double bessel_argument = 2.0 / (3.0 * chi);
  const double bessel_k = std::cyl_bessel_k(1.0 / 3.0, bessel_argument);
  const double energy_ratio =
      constants::electron_rest_energy.eV / photon_energy.eV;

  return Rate{prefactor * energy_ratio * bessel_k * bessel_k};
}

}  // namespace preshaires
