#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "preshaires/conversion_probability.hpp"
#include "preshaires/conversion_sampling.hpp"
#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
#include "preshaires/pair_production.hpp"
#include "preshaires/random.hpp"
#include "preshaires/units.hpp"

namespace {

struct Vector3 {
  double x;
  double y;
  double z;
};

struct ProbabilityInputs {
  preshaires::Energy energy;
  preshaires::StraightTrajectory trajectory;
  std::unique_ptr<preshaires::MagneticFieldModel> field;
};

class Mt19937RandomEngine final : public preshaires::RandomEngine {
 public:
  explicit Mt19937RandomEngine(std::uint64_t seed) : engine_(seed) {}

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

void print_usage() {
  std::cerr
      << "usage:\n"
      << "  preshaires rate --energy-eV VALUE --Bperp-T VALUE\n"
      << "  preshaires probability --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --uniform-field-T Bx,By,Bz\n"
      << "  preshaires probability --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --field-table PATH\n"
      << "  preshaires sample-conversion --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --uniform-field-T Bx,By,Bz --seed VALUE\n"
      << "  preshaires sample-conversion --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --field-table PATH --seed VALUE\n";
}

bool parse_double(std::string_view text, double& value) {
  std::string copy{text};
  char* end = nullptr;
  value = std::strtod(copy.c_str(), &end);
  return end != copy.c_str() && *end == '\0' && std::isfinite(value);
}

Vector3 parse_vector3(std::string_view text) {
  std::stringstream stream{std::string{text}};
  std::string part;
  double values[3]{};
  for (double& value : values) {
    if (!std::getline(stream, part, ',')) {
      throw std::invalid_argument("Expected vector with three comma-separated values");
    }
    if (!parse_double(part, value)) {
      throw std::invalid_argument("Vector contains a nonfinite number");
    }
  }
  if (std::getline(stream, part, ',')) {
    throw std::invalid_argument("Expected vector with exactly three components");
  }
  return Vector3{values[0], values[1], values[2]};
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> columns;
  std::stringstream stream{line};
  std::string column;
  while (std::getline(stream, column, ',')) {
    columns.push_back(column);
  }
  return columns;
}

std::vector<preshaires::TabulatedMagneticFieldNode> read_field_table(
    const std::string& path, preshaires::Length required_length) {
  std::ifstream input{path};
  if (!input) {
    throw std::runtime_error("Cannot open field table");
  }

  std::vector<preshaires::TabulatedMagneticFieldNode> nodes;
  std::string line;
  bool first_line = true;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    if (first_line && line == "s_m,Bx_T,By_T,Bz_T") {
      first_line = false;
      continue;
    }
    first_line = false;

    const auto columns = split_csv_line(line);
    if (columns.size() != 4) {
      throw std::invalid_argument("Field table must have four columns");
    }
    double s = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double bz = 0.0;
    if (!parse_double(columns[0], s) || !parse_double(columns[1], bx) ||
        !parse_double(columns[2], by) || !parse_double(columns[3], bz)) {
      throw std::invalid_argument("Field table contains a nonfinite number");
    }
    nodes.push_back(preshaires::TabulatedMagneticFieldNode{
        preshaires::meters(s),
        preshaires::MagneticFieldVector{preshaires::tesla(bx),
                                        preshaires::tesla(by),
                                        preshaires::tesla(bz)}});
  }

  if (nodes.empty()) {
    throw std::invalid_argument("Field table has no data rows");
  }
  if (nodes.front().path_length.meter > 0.0 ||
      nodes.back().path_length.meter < required_length.meter) {
    throw std::invalid_argument(
        "Field table must cover the complete interval [0, L]");
  }
  return nodes;
}

ProbabilityInputs parse_probability_inputs(int argc, char** argv,
                                           bool allow_seed,
                                           std::uint64_t* seed_out) {
  bool have_energy = false;
  bool have_length = false;
  bool have_direction = false;
  bool have_uniform = false;
  bool have_table = false;
  bool have_seed = false;
  double energy_eV = 0.0;
  double length_m = 0.0;
  Vector3 direction{0.0, 0.0, 0.0};
  Vector3 uniform_field{0.0, 0.0, 0.0};
  std::string table_path;
  std::uint64_t seed = 0;

  for (int i = 2; i < argc; i += 2) {
    if (i + 1 >= argc) {
      throw std::invalid_argument("Missing value for option");
    }
    const std::string_view option{argv[i]};
    if (option == "--energy-eV") {
      have_energy = parse_double(argv[i + 1], energy_eV);
    } else if (option == "--length-m") {
      have_length = parse_double(argv[i + 1], length_m);
    } else if (option == "--direction") {
      direction = parse_vector3(argv[i + 1]);
      have_direction = true;
    } else if (option == "--uniform-field-T") {
      uniform_field = parse_vector3(argv[i + 1]);
      have_uniform = true;
    } else if (option == "--field-table") {
      table_path = argv[i + 1];
      have_table = true;
    } else if (option == "--seed" && allow_seed) {
      double parsed_seed = 0.0;
      have_seed = parse_double(argv[i + 1], parsed_seed) &&
                  parsed_seed >= 0.0 &&
                  parsed_seed <=
                      static_cast<double>(std::numeric_limits<std::uint64_t>::max()) &&
                  std::floor(parsed_seed) == parsed_seed;
      seed = static_cast<std::uint64_t>(parsed_seed);
    } else {
      throw std::invalid_argument("Unknown option");
    }
  }

  if (!have_energy || !have_length || !have_direction ||
      have_uniform == have_table || (allow_seed && !have_seed)) {
    throw std::invalid_argument("Missing or inconsistent probability arguments");
  }

  const auto length = preshaires::meters(length_m);
  auto field = std::unique_ptr<preshaires::MagneticFieldModel>{};
  if (have_uniform) {
    field = std::make_unique<preshaires::UniformMagneticField>(
        preshaires::MagneticFieldVector{preshaires::tesla(uniform_field.x),
                                        preshaires::tesla(uniform_field.y),
                                        preshaires::tesla(uniform_field.z)});
  } else {
    field = std::make_unique<preshaires::TabulatedMagneticField>(
        read_field_table(table_path, length));
  }

  if (seed_out != nullptr) {
    *seed_out = seed;
  }
  return ProbabilityInputs{
      preshaires::electron_volts(energy_eV),
      preshaires::StraightTrajectory{
          preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                               preshaires::meters(0.0)},
          preshaires::make_direction(direction.x, direction.y, direction.z),
          length},
      std::move(field)};
}

int run_rate(int argc, char** argv) {
  if (argc != 6) {
    print_usage();
    return 2;
  }

  bool have_energy = false;
  bool have_bperp = false;
  double energy_eV = 0.0;
  double bperp_T = 0.0;

  for (int i = 2; i < argc; i += 2) {
    const std::string_view option{argv[i]};
    double value = 0.0;
    if (!parse_double(argv[i + 1], value)) {
      print_usage();
      return 2;
    }

    if (option == "--energy-eV") {
      energy_eV = value;
      have_energy = true;
    } else if (option == "--Bperp-T") {
      bperp_T = value;
      have_bperp = true;
    } else {
      print_usage();
      return 2;
    }
  }

  if (!have_energy || !have_bperp) {
    print_usage();
    return 2;
  }

  const auto energy = preshaires::electron_volts(energy_eV);
  const auto field = preshaires::tesla(bperp_T);
  const double chi = preshaires::photon_chi(energy, field);
  const auto rate = preshaires::erber1966_pair_production_rate(energy, field);

  std::cout << std::setprecision(17) << "model Erber1966\n"
            << "chi " << chi << '\n'
            << "rate_per_second " << rate.per_second << '\n';

  return 0;
}

int run_probability(int argc, char** argv) {
  const ProbabilityInputs inputs =
      parse_probability_inputs(argc, argv, false, nullptr);
  const preshaires::IntegrationOptions options{1.0e-8, 1.0e-14, 100000};

  const preshaires::ConversionProbabilityResult result =
      preshaires::conversion_probability(inputs.energy, inputs.trajectory,
                                         *inputs.field, options);

  std::cout << std::setprecision(17) << "model Erber1966\n"
            << "optical_depth " << result.optical_depth.value << '\n'
            << "probability " << result.probability.value << '\n'
            << "evaluations " << result.evaluations << '\n';

  return 0;
}

int run_sample_conversion(int argc, char** argv) {
  std::uint64_t seed = 0;
  const ProbabilityInputs inputs =
      parse_probability_inputs(argc, argv, true, &seed);
  Mt19937RandomEngine rng{seed};
  const preshaires::IntegrationOptions integration_options{1.0e-8, 1.0e-14,
                                                          100000};
  const preshaires::LocalizationOptions localization_options{
      preshaires::meters(1.0e-3), 1.0e-12, 128};
  const preshaires::ConversionSample sample =
      preshaires::sample_first_conversion(inputs.energy, inputs.trajectory,
                                          *inputs.field, rng,
                                          integration_options,
                                          localization_options);

  std::cout << std::setprecision(17) << "model Erber1966\n"
            << "converted " << (sample.converted ? "true" : "false") << '\n'
            << "random_uniform " << sample.random_uniform << '\n'
            << "target_optical_depth " << sample.target_optical_depth.value
            << '\n'
            << "total_optical_depth " << sample.total_optical_depth.value
            << '\n'
            << "total_probability " << sample.total_probability.value << '\n'
            << "evaluations " << sample.evaluations << '\n';
  if (sample.converted) {
    const auto position = *sample.position;
    std::cout << "path_length_m " << sample.path_length->meter << '\n'
              << "position_m " << position.x.meter << ',' << position.y.meter
              << ',' << position.z.meter << '\n';
  }

  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      print_usage();
      return 2;
    }
    const std::string_view command{argv[1]};
    if (command == "rate") {
      return run_rate(argc, argv);
    }
    if (command == "probability") {
      return run_probability(argc, argv);
    }
    if (command == "sample-conversion") {
      return run_sample_conversion(argc, argv);
    }
    print_usage();
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
