#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>

#include "preshaires/pair_production.hpp"
#include "preshaires/units.hpp"

namespace {

void print_usage() {
  std::cerr << "usage: preshaires rate --energy-eV VALUE --Bperp-T VALUE\n";
}

bool parse_double(std::string_view text, double& value) {
  std::string copy{text};
  char* end = nullptr;
  value = std::strtod(copy.c_str(), &end);
  return end != copy.c_str() && *end == '\0';
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 6 || std::string_view{argv[1]} != "rate") {
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
