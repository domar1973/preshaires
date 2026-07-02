#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "preshaires/conversion_probability.hpp"
#include "preshaires/geometry.hpp"
#include "preshaires/magnetic_field.hpp"
#include "preshaires/pair_production.hpp"
#include "preshaires/units.hpp"

namespace {

struct Vector3 {
  double x;
  double y;
  double z;
};

void print_usage() {
  std::cerr
      << "usage:\n"
      << "  preshaires rate --energy-eV VALUE --Bperp-T VALUE\n"
      << "  preshaires probability --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --uniform-field-T Bx,By,Bz\n"
      << "  preshaires probability --energy-eV VALUE --length-m VALUE "
         "--direction X,Y,Z --field-table PATH\n";
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
  bool have_energy = false;
  bool have_length = false;
  bool have_direction = false;
  bool have_uniform = false;
  bool have_table = false;
  double energy_eV = 0.0;
  double length_m = 0.0;
  Vector3 direction{0.0, 0.0, 0.0};
  Vector3 uniform_field{0.0, 0.0, 0.0};
  std::string table_path;

  for (int i = 2; i < argc; i += 2) {
    if (i + 1 >= argc) {
      print_usage();
      return 2;
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
    } else {
      print_usage();
      return 2;
    }
  }

  if (!have_energy || !have_length || !have_direction ||
      have_uniform == have_table) {
    print_usage();
    return 2;
  }

  const auto length = preshaires::meters(length_m);
  const preshaires::StraightTrajectory trajectory{
      preshaires::Position{preshaires::meters(0.0), preshaires::meters(0.0),
                           preshaires::meters(0.0)},
      preshaires::make_direction(direction.x, direction.y, direction.z), length};
  const auto energy = preshaires::electron_volts(energy_eV);
  const preshaires::IntegrationOptions options{1.0e-8, 1.0e-14, 100000};

  preshaires::ConversionProbabilityResult result{};
  if (have_uniform) {
    const preshaires::UniformMagneticField field{
        preshaires::MagneticFieldVector{preshaires::tesla(uniform_field.x),
                                        preshaires::tesla(uniform_field.y),
                                        preshaires::tesla(uniform_field.z)}};
    result =
        preshaires::conversion_probability(energy, trajectory, field, options);
  } else {
    auto nodes = read_field_table(table_path, length);
    const preshaires::TabulatedMagneticField field{std::move(nodes)};
    result =
        preshaires::conversion_probability(energy, trajectory, field, options);
  }

  std::cout << std::setprecision(17) << "model Erber1966\n"
            << "optical_depth " << result.optical_depth.value << '\n'
            << "probability " << result.probability.value << '\n'
            << "evaluations " << result.evaluations << '\n';

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
    print_usage();
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
