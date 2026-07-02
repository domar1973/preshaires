#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "preshaires/conversion_probability.hpp"
#include "preshaires/conversion_sampling.hpp"
#include "preshaires/constants.hpp"
#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
#include "preshaires/pair_production.hpp"
#include "preshaires/random.hpp"
#include "preshaires/units.hpp"

namespace {

int failures = 0;

class FixedRandomEngine final : public preshaires::RandomEngine {
 public:
  explicit FixedRandomEngine(double value) : value_(value) {}
  double uniform_open01() override { return value_; }

 private:
  double value_;
};

class SequenceRandomEngine final : public preshaires::RandomEngine {
 public:
  explicit SequenceRandomEngine(std::uint64_t seed) : engine_(seed) {}
  double uniform_open01() override {
    for (;;) {
      const std::uint64_t value = engine_();
      if (value != 0 && value != std::numeric_limits<std::uint64_t>::max()) {
        return static_cast<double>(value) /
               static_cast<double>(std::numeric_limits<std::uint64_t>::max());
      }
    }
  }

 private:
  std::mt19937_64 engine_;
};

class StepMagneticField final : public preshaires::MagneticFieldModel {
 public:
  StepMagneticField(preshaires::Length switch_path,
                    preshaires::MagneticFieldVector before,
                    preshaires::MagneticFieldVector after)
      : switch_path_(switch_path), before_(before), after_(after) {}

  preshaires::MagneticFieldVector field_at(const preshaires::Position&,
                                           preshaires::Length path) const override {
    return path.meter < switch_path_.meter ? before_ : after_;
  }

 private:
  preshaires::Length switch_path_;
  preshaires::MagneticFieldVector before_;
  preshaires::MagneticFieldVector after_;
};

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

  const auto tau = preshaires::optical_depth(energy, trajectory, field, length,
                                             options);
  expect_near(tau.optical_depth.value, result.optical_depth.value, 0.0,
              "optical_depth preserves v0.2 tau");
  expect_near(preshaires::conversion_probability_from_optical_depth(
                  tau.optical_depth)
                  .value,
              result.probability.value, 0.0,
              "conversion_probability preserves v0.2 probability");
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

preshaires::IntegrationOptions sampling_integration_options() {
  return preshaires::IntegrationOptions{1.0e-10, 1.0e-14, 200000};
}

preshaires::LocalizationOptions sampling_localization_options() {
  return preshaires::LocalizationOptions{preshaires::meters(1.0e-5), 1.0e-12,
                                         160};
}

void test_uniform_sampling_conversion() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto length = preshaires::meters(1.0e7);
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), length};
  const auto b_perp = preshaires::tesla(2.115138e-5);
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{b_perp, preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  constexpr double uniform = 0.5;
  FixedRandomEngine rng{uniform};
  const auto sample = preshaires::sample_first_conversion(
      energy, trajectory, field, rng, sampling_integration_options(),
      sampling_localization_options());
  const double alpha =
      preshaires::erber1966_pair_production_rate(energy, b_perp).per_second /
      preshaires::constants::speed_of_light_m_per_s;
  const double expected_s = -std::log(uniform) / alpha;
  expect_true(sample.converted, "uniform sample converts");
  expect_near(sample.path_length->meter, expected_s, 1.0e-10,
              "uniform sample path length");
  const auto expected_position =
      preshaires::position_at(trajectory, preshaires::meters(expected_s));
  expect_near(sample.position->z.meter, expected_position.z.meter, 1.0e-10,
              "uniform sample position");
  const auto tau_at_sample = preshaires::optical_depth(
      energy, trajectory, field, *sample.path_length,
      sampling_integration_options());
  expect_near(tau_at_sample.optical_depth.value,
              sample.target_optical_depth.value, 1.0e-8,
              "tau at sampled path matches target");
}

void test_sampling_no_conversion_and_parallel() {
  const auto trajectory = make_test_trajectory(1.0e7);
  const auto energy = preshaires::electron_volts(7.0e19);
  const preshaires::UniformMagneticField perpendicular{
      preshaires::MagneticFieldVector{preshaires::tesla(2.115138e-5),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  FixedRandomEngine no_conversion_rng{1.0e-6};
  const auto no_conversion = preshaires::sample_first_conversion(
      energy, trajectory, perpendicular, no_conversion_rng,
      sampling_integration_options(), sampling_localization_options());
  expect_true(!no_conversion.converted, "large target optical depth does not convert");
  expect_true(!no_conversion.path_length.has_value(), "no conversion has no path");
  expect_true(!no_conversion.position.has_value(), "no conversion has no position");

  const preshaires::UniformMagneticField parallel{
      preshaires::MagneticFieldVector{preshaires::tesla(0.0),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(2.0e-5)}};
  FixedRandomEngine parallel_rng{0.5};
  const auto parallel_sample = preshaires::sample_first_conversion(
      energy, trajectory, parallel, parallel_rng, sampling_integration_options(),
      sampling_localization_options());
  expect_true(!parallel_sample.converted, "parallel field never converts");
  expect_near(parallel_sample.total_optical_depth.value, 0.0, 0.0,
              "parallel total tau is zero");
}

void test_sampling_near_extremes() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto length = preshaires::meters(1.0e7);
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), length};
  const auto b_perp = preshaires::tesla(2.115138e-5);
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{b_perp, preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  FixedRandomEngine near_start_rng{std::nextafter(1.0, 0.0)};
  const auto near_start = preshaires::sample_first_conversion(
      energy, trajectory, field, near_start_rng, sampling_integration_options(),
      sampling_localization_options());
  expect_true(near_start.converted, "near-start sample converts");
  expect_true(near_start.path_length->meter < 1.0,
              "near-start sample is close to start");

  const auto total = preshaires::optical_depth(
      energy, trajectory, field, length, sampling_integration_options());
  const double target = std::nextafter(total.optical_depth.value, 0.0);
  FixedRandomEngine near_end_rng{std::exp(-target)};
  const auto near_end = preshaires::sample_first_conversion(
      energy, trajectory, field, near_end_rng, sampling_integration_options(),
      sampling_localization_options());
  expect_true(near_end.converted, "near-end sample converts");
  expect_true(near_end.path_length->meter <= length.meter,
              "near-end sample stays inside trajectory");
  expect_true(near_end.path_length->meter > length.meter - 100.0,
              "near-end sample is close to end");
}

void test_invalid_rng_values() {
  const auto trajectory = make_test_trajectory(1.0e7);
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{preshaires::tesla(2.0e-5),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  const auto call_with = [&](double value) {
    FixedRandomEngine rng{value};
    (void)preshaires::sample_first_conversion(
        preshaires::electron_volts(7.0e19), trajectory, field, rng,
        sampling_integration_options(), sampling_localization_options());
  };
  expect_throws<std::runtime_error>([&] { call_with(0.0); }, "rng zero rejected");
  expect_throws<std::runtime_error>([&] { call_with(1.0); }, "rng one rejected");
  expect_throws<std::runtime_error>([&] { call_with(-0.1); }, "rng negative rejected");
  expect_throws<std::runtime_error>(
      [&] { call_with(std::numeric_limits<double>::quiet_NaN()); },
      "rng NaN rejected");
  expect_throws<std::runtime_error>(
      [&] { call_with(std::numeric_limits<double>::infinity()); },
      "rng infinity rejected");
}

void test_sampling_zero_rate_interval_and_plateau() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), preshaires::meters(1000.0)};
  const auto active_field = preshaires::MagneticFieldVector{
      preshaires::tesla(1.0e-4), preshaires::tesla(0.0),
      preshaires::tesla(0.0)};
  const StepMagneticField delayed{
      preshaires::meters(200.0),
      preshaires::MagneticFieldVector{preshaires::tesla(0.0),
                                      preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)},
      active_field};
  const double alpha = preshaires::erber1966_pair_production_rate(
                           energy, preshaires::tesla(1.0e-4))
                           .per_second /
                       preshaires::constants::speed_of_light_m_per_s;
  const double target = alpha * 50.0;
  FixedRandomEngine delayed_rng{std::exp(-target)};
  const auto delayed_sample = preshaires::sample_first_conversion(
      energy, trajectory, delayed, delayed_rng,
      preshaires::IntegrationOptions{1.0e-9, 1.0e-13, 200000},
      sampling_localization_options());
  expect_true(delayed_sample.converted, "delayed field sample converts");
  expect_true(delayed_sample.path_length->meter >= 200.0,
              "sample is not in zero-rate interval");

  const preshaires::TabulatedMagneticField plateau{
      std::vector<preshaires::TabulatedMagneticFieldNode>{
          {preshaires::meters(0.0), active_field},
          {preshaires::meters(500.0),
           {preshaires::tesla(0.0), preshaires::tesla(0.0),
            preshaires::tesla(0.0)}},
          {preshaires::meters(700.0),
           {preshaires::tesla(0.0), preshaires::tesla(0.0),
            preshaires::tesla(0.0)}},
          {preshaires::meters(1000.0), active_field}}};
  const double plateau_target = preshaires::optical_depth(
                                    energy, trajectory, plateau,
                                    preshaires::meters(500.0),
                                    preshaires::IntegrationOptions{
                                        1.0e-9, 1.0e-13, 200000})
                                    .optical_depth.value;
  FixedRandomEngine plateau_rng{std::exp(-plateau_target)};
  const auto plateau_sample = preshaires::sample_first_conversion(
      energy, trajectory, plateau, plateau_rng,
      preshaires::IntegrationOptions{1.0e-9, 1.0e-13, 200000},
      preshaires::LocalizationOptions{preshaires::meters(2.0), 0.0, 160});
  expect_true(plateau_sample.converted, "plateau sample converts");
  expect_near(plateau_sample.path_length->meter, 500.0, 5.0e-3,
              "plateau sample returns first compatible point");
}

void test_sampling_statistics_and_reproducibility() {
  const auto energy = preshaires::electron_volts(7.0e19);
  const auto length = preshaires::meters(1.0e6);
  const auto trajectory = preshaires::StraightTrajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(0.0, 0.0, -1.0), length};
  const auto b_perp = preshaires::tesla(2.115138e-5);
  const preshaires::UniformMagneticField field{
      preshaires::MagneticFieldVector{b_perp, preshaires::tesla(0.0),
                                      preshaires::tesla(0.0)}};
  const double alpha =
      preshaires::erber1966_pair_production_rate(energy, b_perp).per_second /
      preshaires::constants::speed_of_light_m_per_s;
  const double expected_fraction = -std::expm1(-alpha * length.meter);
  constexpr int samples = 10000;
  int converted = 0;
  SequenceRandomEngine rng{123456};
  for (int i = 0; i < samples; ++i) {
    const auto sample = preshaires::sample_first_conversion(
        energy, trajectory, field, rng,
        preshaires::IntegrationOptions{1.0e-8, 1.0e-12, 10000},
        preshaires::LocalizationOptions{preshaires::meters(1.0e-2), 1.0e-10,
                                        80});
    if (sample.converted) {
      ++converted;
    }
  }
  const double observed_fraction = static_cast<double>(converted) / samples;
  const double sigma =
      std::sqrt(expected_fraction * (1.0 - expected_fraction) / samples);
  expect_true(std::abs(observed_fraction - expected_fraction) < 5.0 * sigma,
              "sample conversion fraction matches binomial expectation");

  SequenceRandomEngine rng_a{98765};
  SequenceRandomEngine rng_b{98765};
  for (int i = 0; i < 10; ++i) {
    const auto a = preshaires::sample_first_conversion(
        energy, trajectory, field, rng_a,
        preshaires::IntegrationOptions{1.0e-8, 1.0e-12, 10000},
        sampling_localization_options());
    const auto b = preshaires::sample_first_conversion(
        energy, trajectory, field, rng_b,
        preshaires::IntegrationOptions{1.0e-8, 1.0e-12, 10000},
        sampling_localization_options());
    expect_true(a.converted == b.converted, "reproducible converted flag");
    expect_near(a.random_uniform, b.random_uniform, 0.0,
                "reproducible random U");
    expect_near(a.target_optical_depth.value, b.target_optical_depth.value, 0.0,
                "reproducible target tau");
    if (a.converted) {
      expect_near(a.path_length->meter, b.path_length->meter, 0.0,
                  "reproducible path length");
    }
  }
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
  test_uniform_sampling_conversion();
  test_sampling_no_conversion_and_parallel();
  test_sampling_near_extremes();
  test_invalid_rng_values();
  test_sampling_zero_rate_interval_and_plateau();
  test_sampling_statistics_and_reproducibility();

  if (failures != 0) {
    std::cerr << failures << " test failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
