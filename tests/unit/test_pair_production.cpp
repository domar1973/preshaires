#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "preshaires/constants.hpp"
#include "preshaires/pair_production.hpp"
#include "preshaires/units.hpp"

namespace {

int failures = 0;

void expect_true(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void expect_near(double actual, double expected, double rel_tol,
                 const std::string& message) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  if (std::abs(actual - expected) > rel_tol * scale) {
    std::cerr << "FAIL: " << message << ": actual=" << actual
              << " expected=" << expected << '\n';
    ++failures;
  }
}

double reference_erber_rate(preshaires::Energy energy,
                            preshaires::MagneticField field) {
  constexpr double coefficient = 0.16;
  const double energy_joule =
      energy.eV * preshaires::constants::joule_per_electron_volt;
  const double electron_rest_energy_joule =
      preshaires::constants::electron_rest_energy.eV *
      preshaires::constants::joule_per_electron_volt;
  const double chi = preshaires::photon_chi(energy, field);
  if (chi <= 0.0) {
    return 0.0;
  }
  const double prefactor =
      coefficient * preshaires::constants::fine_structure_constant *
      electron_rest_energy_joule /
      preshaires::constants::reduced_planck_constant_J_s;
  const double bessel_k = std::cyl_bessel_k(1.0 / 3.0, 2.0 / (3.0 * chi));
  return prefactor * (electron_rest_energy_joule / energy_joule) * bessel_k *
         bessel_k;
}

void test_zero_transverse_field() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto field = preshaires::tesla(0.0);
  expect_near(preshaires::photon_chi(energy, field), 0.0, 0.0,
              "B_perp = 0 gives chi = 0");
  expect_near(
      preshaires::erber1966_pair_production_rate(energy, field).per_second,
      0.0, 0.0, "B_perp = 0 gives rate = 0");
}

void test_diagnostic_case() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto field = preshaires::tesla(2.115138e-5);
  const double chi = preshaires::photon_chi(energy, field);
  expect_near(chi, 0.3281999188925183, 2.0e-10,
              "diagnostic chi matches reference");
  expect_near(
      preshaires::erber1966_pair_production_rate(energy, field).per_second,
      reference_erber_rate(energy, field), 1.0e-14,
      "diagnostic rate matches std::cyl_bessel_k reference");
}

void test_positivity_monotonicity_and_finiteness() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto low_field = preshaires::tesla(2.0e-5);
  const auto high_field = preshaires::tesla(2.1e-5);
  const auto low_rate =
      preshaires::erber1966_pair_production_rate(energy, low_field);
  const auto high_rate =
      preshaires::erber1966_pair_production_rate(energy, high_field);
  expect_true(low_rate.per_second > 0.0, "rate is positive");
  expect_true(high_rate.per_second > low_rate.per_second,
              "rate is locally monotonic in B_perp");
  expect_true(std::isfinite(low_rate.per_second), "rate is finite");
  expect_true(std::isfinite(high_rate.per_second), "higher rate is finite");
}

void test_no_nan_on_grid() {
  for (double energy_eV : {1.0e18, 7.0e19, 1.0e21}) {
    for (double field_T : {0.0, 1.0e-9, 2.0e-5, 1.0e-4}) {
      const auto rate = preshaires::erber1966_pair_production_rate(
          preshaires::electron_volts(energy_eV), preshaires::tesla(field_T));
      expect_true(!std::isnan(rate.per_second), "grid rate is not NaN");
    }
  }
}

void test_unit_conversions() {
  expect_near(preshaires::gigaelectron_volts(70.0).eV, 7.0e10, 0.0,
              "GeV to eV conversion");
  expect_near(preshaires::nanotesla(21151.38).tesla, 2.115138e-5, 1.0e-15,
              "nT to T conversion");
}

void test_metadata() {
  const auto metadata = preshaires::erber1966_metadata();
  expect_true(metadata.model == preshaires::PairProductionModel::Erber1966,
              "metadata model");
  expect_true(!metadata.reference.empty(), "metadata reference");
  expect_true(!metadata.equations.empty(), "metadata equations");
  expect_true(!metadata.assumptions.empty(), "metadata assumptions");
  expect_true(!metadata.validity_range.empty(), "metadata validity range");
}

}  // namespace

int main() {
  test_zero_transverse_field();
  test_diagnostic_case();
  test_positivity_monotonicity_and_finiteness();
  test_no_nan_on_grid();
  test_unit_conversions();
  test_metadata();

  if (failures != 0) {
    std::cerr << failures << " test failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
