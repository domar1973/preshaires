#pragma once

#include "preshaires/units.hpp"

namespace preshaires::constants {

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double speed_of_light_m_per_s = 299792458.0;
inline constexpr double reduced_planck_constant_J_s = 1.05457e-34;
inline constexpr double elementary_charge_C = 1.6021892e-19;
inline constexpr double joule_per_electron_volt = elementary_charge_C;
inline constexpr double electron_mass_kg = 9.109534e-31;
inline constexpr double fine_structure_constant = 1.0 / 137.0;

inline constexpr Energy electron_rest_energy{
    electron_mass_kg * speed_of_light_m_per_s * speed_of_light_m_per_s /
    joule_per_electron_volt};

inline constexpr MagneticField critical_magnetic_field{
    electron_mass_kg * electron_mass_kg * speed_of_light_m_per_s *
    speed_of_light_m_per_s /
    (elementary_charge_C * reduced_planck_constant_J_s)};

}  // namespace preshaires::constants
