#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "preshaires/conversion_probability.hpp"
#include "preshaires/constants.hpp"
#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
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

template <typename Exception, typename Callable>
void expect_throws(Callable callable, const std::string& message) {
  try {
    callable();
  } catch (const Exception&) {
    return;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << message << ": wrong exception: " << error.what()
              << '\n';
    ++failures;
    return;
  }
  std::cerr << "FAIL: " << message << ": no exception\n";
  ++failures;
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

preshaires::StraightTrajectory make_test_trajectory(double length_m) {
  return preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(1.0), preshaires::meters(2.0),
                           preshaires::meters(3.0)},
      preshaires::make_direction(0.0, 0.0, -1.0),
      preshaires::meters(length_m)};
}

void test_geometry() {
  const auto direction = preshaires::make_direction(3.0, 4.0, 0.0);
  expect_near(std::hypot(direction.x, direction.y, direction.z), 1.0, 1.0e-15,
              "direction is normalized");
  expect_throws<std::invalid_argument>(
      [] { (void)preshaires::make_direction(0.0, 0.0, 0.0); },
      "zero direction is rejected");

  const auto trajectory = make_test_trajectory(10.0);
  const auto position = preshaires::position_at(trajectory, preshaires::meters(4.0));
  expect_near(position.x.meter, 1.0, 0.0, "position x on trajectory");
  expect_near(position.y.meter, 2.0, 0.0, "position y on trajectory");
  expect_near(position.z.meter, -1.0, 0.0, "position z on trajectory");

  const auto along_z = preshaires::make_direction(0.0, 0.0, 1.0);
  const auto parallel = preshaires::transverse_field(
      preshaires::MagneticFieldVector{preshaires::tesla(0.0),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(2.0e-5)},
      along_z);
  expect_near(parallel.tesla, 0.0, 0.0, "parallel B gives zero B_perp");

  const auto perpendicular = preshaires::transverse_field(
      preshaires::MagneticFieldVector{preshaires::tesla(2.0e-5),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)},
      along_z);
  expect_near(perpendicular.tesla, 2.0e-5, 1.0e-15,
              "perpendicular B gives magnitude");

  const auto oblique = preshaires::transverse_field(
      preshaires::MagneticFieldVector{preshaires::tesla(3.0),
                                      preshaires::tesla(4.0),
                                      preshaires::tesla(0.0)},
      preshaires::make_direction(1.0, 0.0, 0.0));
  expect_near(oblique.tesla, 4.0, 1.0e-15, "oblique B_perp");
}

void test_uniform_field() {
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{preshaires::tesla(1.0e-5),
                                      preshaires::tesla(2.0e-5),
                                      preshaires::tesla(3.0e-5)}};
  const auto first = field.field_at(
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::meters(0.0));
  const auto second = field.field_at(
      preshaires::Position{preshaires::meters(9.0), preshaires::meters(8.0),
                           preshaires::meters(7.0)},
      preshaires::meters(42.0));
  expect_near(first.x.tesla, second.x.tesla, 0.0, "uniform field x");
  expect_near(first.y.tesla, second.y.tesla, 0.0, "uniform field y");
  expect_near(first.z.tesla, second.z.tesla, 0.0, "uniform field z");
}

std::vector<preshaires::TabulatedMagneticFieldNode> tabulated_nodes() {
  return {
      {preshaires::meters(0.0),
       {preshaires::tesla(0.0), preshaires::tesla(1.0e-5),
        preshaires::tesla(0.0)}},
      {preshaires::meters(10.0),
       {preshaires::tesla(2.0e-5), preshaires::tesla(3.0e-5),
        preshaires::tesla(4.0e-5)}}};
}

void test_tabulated_field() {
  const preshaires::Position origin{preshaires::meters(0.0),
                                    preshaires::meters(0.0),
                                    preshaires::meters(0.0)};
  const preshaires::TabulatedMagneticField field{tabulated_nodes()};
  const auto left = field.field_at(origin, preshaires::meters(0.0));
  const auto middle = field.field_at(origin, preshaires::meters(5.0));
  const auto right = field.field_at(origin, preshaires::meters(10.0));
  expect_near(left.y.tesla, 1.0e-5, 0.0, "tabulated left node exact");
  expect_near(right.z.tesla, 4.0e-5, 0.0, "tabulated right node exact");
  expect_near(middle.x.tesla, 1.0e-5, 1.0e-15, "tabulated interpolation x");
  expect_near(middle.y.tesla, 2.0e-5, 1.0e-15, "tabulated interpolation y");
  expect_near(middle.z.tesla, 2.0e-5, 1.0e-15, "tabulated interpolation z");

  expect_throws<std::invalid_argument>(
      [] {
        (void)preshaires::TabulatedMagneticField{
            std::vector<preshaires::TabulatedMagneticFieldNode>{
                tabulated_nodes()[0], tabulated_nodes()[0]}};
      },
      "tabulated repeated nodes rejected");
  expect_throws<std::invalid_argument>(
      [] {
        auto nodes = tabulated_nodes();
        std::swap(nodes[0], nodes[1]);
        (void)preshaires::TabulatedMagneticField{nodes};
      },
      "tabulated incorrect order rejected");
  expect_throws<std::out_of_range>(
      [&] { (void)field.field_at(origin, preshaires::meters(11.0)); },
      "tabulated extrapolation rejected");
  expect_throws<std::invalid_argument>(
      [] {
        auto nodes = tabulated_nodes();
        nodes[1].field.x = preshaires::tesla(
            std::numeric_limits<double>::quiet_NaN());
        (void)preshaires::TabulatedMagneticField{nodes};
      },
      "tabulated NaN rejected");
}

void test_uniform_conversion_probability() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto length = preshaires::meters(1.0e7);
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), length};
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{preshaires::tesla(2.115138e-5),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  const preshaires::IntegrationOptions options{1.0e-10, 1.0e-14, 10000};
  const auto result =
      preshaires::conversion_probability(energy, trajectory, field, options);
  const double rate = preshaires::erber1966_pair_production_rate(
                          energy, preshaires::tesla(2.115138e-5))
                          .per_second;
  const double tau_expected =
      rate * length.meter / preshaires::constants::speed_of_light_m_per_s;
  expect_near(result.optical_depth.value, tau_expected, 1.0e-13,
              "uniform optical depth analytic");
  expect_near(result.probability.value, -std::expm1(-tau_expected), 1.0e-13,
              "uniform probability analytic");
}

void test_parallel_conversion_probability() {
  const auto trajectory = make_test_trajectory(1.0e7);
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{preshaires::tesla(0.0),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(2.0e-5)}};
  const auto result = preshaires::conversion_probability(
      preshaires::electron_volts(7.0e19), trajectory, field,
      preshaires::IntegrationOptions{1.0e-10, 1.0e-14, 10000});
  expect_near(result.optical_depth.value, 0.0, 0.0, "parallel tau zero");
  expect_near(result.probability.value, 0.0, 0.0, "parallel P zero");
}

double high_resolution_reference(preshaires::Energy energy, double length_m) {
  constexpr int intervals = 20000;
  const double h = length_m / intervals;
  double sum = 0.0;
  for (int i = 0; i <= intervals; ++i) {
    const double s = i * h;
    const double b = (1.0e-5 + 1.0e-5 * s / length_m);
    const double value =
        preshaires::erber1966_pair_production_rate(energy, preshaires::tesla(b))
            .per_second /
        preshaires::constants::speed_of_light_m_per_s;
    const double weight = (i == 0 || i == intervals) ? 0.5 : 1.0;
    sum += weight * value;
  }
  return sum * h;
}

void test_linear_profile_and_convergence() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const double length_m = 1.0e7;
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), preshaires::meters(length_m)};
  const preshaires::TabulatedMagneticField field{
      std::vector<preshaires::TabulatedMagneticFieldNode>{
          {preshaires::meters(0.0),
           {preshaires::tesla(1.0e-5), preshaires::tesla(0.0),
            preshaires::tesla(0.0)}},
          {preshaires::meters(length_m),
           {preshaires::tesla(2.0e-5), preshaires::tesla(0.0),
            preshaires::tesla(0.0)}}}};
  const auto loose = preshaires::conversion_probability(
      energy, trajectory, field,
      preshaires::IntegrationOptions{1.0e-5, 1.0e-10, 100000});
  const auto tight = preshaires::conversion_probability(
      energy, trajectory, field,
      preshaires::IntegrationOptions{1.0e-9, 1.0e-14, 100000});
  const double reference = high_resolution_reference(energy, length_m);
  expect_near(tight.optical_depth.value, reference, 2.0e-7,
              "linear profile high-resolution reference");
  expect_near(tight.optical_depth.value, loose.optical_depth.value, 5.0e-4,
              "tight tolerance stabilizes result");
}

void test_small_and_saturated_probability() {
  const auto tiny_tau_result = preshaires::ConversionProbabilityResult{
      preshaires::OpticalDepth{0.0}, preshaires::Probability{0.0}, 0};
  expect_near(tiny_tau_result.probability.value, 0.0, 0.0,
              "explicit zero probability fixture");

  const auto energy = preshaires::electron_volts(7.0e19);
  const auto direction = preshaires::make_direction(0.0, 0.0, -1.0);
  const auto small_trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      direction, preshaires::meters(1.0e-6)};
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{preshaires::tesla(2.115138e-5),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  const auto small = preshaires::conversion_probability(
      energy, small_trajectory, field,
      preshaires::IntegrationOptions{1.0e-12, 1.0e-25, 1000});
  expect_near(small.probability.value, small.optical_depth.value, 1.0e-12,
              "small probability preserves tau precision");

  const auto long_trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      direction, preshaires::meters(1.0e10)};
  const preshaires::UniformMagneticField strong_field{
      preshaires::MagneticFieldVector{preshaires::tesla(1.0e-3),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  const auto saturated = preshaires::conversion_probability(
      energy, long_trajectory, strong_field,
      preshaires::IntegrationOptions{1.0e-8, 1.0e-12, 10000});
  expect_true(saturated.probability.value <= 1.0,
              "saturated probability does not exceed one");
  expect_true(saturated.probability.value > 0.999999,
              "saturated probability approaches one");
}

}  // namespace

int main() {
  test_zero_transverse_field();
  test_diagnostic_case();
  test_positivity_monotonicity_and_finiteness();
  test_no_nan_on_grid();
  test_unit_conversions();
  test_metadata();
  test_geometry();
  test_uniform_field();
  test_tabulated_field();
  test_uniform_conversion_probability();
  test_parallel_conversion_probability();
  test_linear_profile_and_convergence();
  test_small_and_saturated_probability();

  if (failures != 0) {
    std::cerr << failures << " test failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
