# preshaires Handoff v0.3

This note summarizes the current development state so work can resume in a new
session without reconstructing context from the full history.

## 1. Estado Alcanzado

- v0.1: local Erber 1966 pair-production rate for a photon in a transverse
  magnetic field.
- v0.2: optical depth and deterministic conversion probability over a
  prescribed straight trajectory.
- v0.3: first-conversion sampling by optical depth, using
  `tau_target = -log(U)`.

## 2. Arquitectura Actual

CMake targets:

- `preshaires-core`: C++20 physical core.
- `preshaires-cli`: standalone command-line interface.
- `preshaires-aires-adapter`: placeholder adapter, no real AIRES linkage yet.
- `preshaires-tests`: unit and lightweight CLI regression tests.

Main public headers:

- `include/preshaires/units.hpp`
- `include/preshaires/constants.hpp`
- `include/preshaires/pair_production.hpp`
- `include/preshaires/geometry.hpp`
- `include/preshaires/magnetic_field.hpp`
- `include/preshaires/conversion_probability.hpp`
- `include/preshaires/random.hpp`
- `include/preshaires/conversion_sampling.hpp`

Separation:

- The core owns physics, unit wrappers, geometry primitives, field interfaces,
  optical-depth integration and conversion sampling.
- The CLI owns argument parsing, CSV field-table reading, RNG construction and
  output formatting.
- AIRES must remain confined to `adapters/aires`; no AIRES headers or symbols
  are used in the core.

`RandomEngine` contract:

```cpp
class RandomEngine {
 public:
  virtual ~RandomEngine() = default;
  virtual double uniform_open01() = 0;
};
```

`uniform_open01()` must return a finite value `0 < U < 1`. The core rejects
`NaN`, infinities, `U <= 0` and `U >= 1` explicitly. There is no global RNG.

The core currently has no mutable global state and no external dependencies
beyond the C++20 standard library.

## 3. API Publica Actual

Important unit/value types:

```cpp
struct Energy { double eV; };
struct MagneticField { double tesla; };
struct Length { double meter; };
struct Time { double second; };
struct Rate { double per_second; };

struct OpticalDepth { double value; };
struct Probability { double value; };
```

Geometry and fields:

```cpp
struct Position;
struct Direction;
struct MagneticFieldVector;
struct StraightTrajectory;

Direction make_direction(double x, double y, double z);
Position position_at(const StraightTrajectory& trajectory, Length path_length);
MagneticField transverse_field(const MagneticFieldVector& field,
                               const Direction& direction);

class MagneticFieldModel {
 public:
  virtual ~MagneticFieldModel() = default;
  virtual MagneticFieldVector field_at(const Position& position,
                                       Length path_length) const = 0;
};

class UniformMagneticField;
class TabulatedMagneticField;
```

Local rate:

```cpp
double photon_chi(Energy photon_energy,
                  MagneticField transverse_magnetic_field);

Rate erber1966_pair_production_rate(
    Energy photon_energy,
    MagneticField transverse_magnetic_field);
```

Integration:

```cpp
struct IntegrationOptions {
  double relative_tolerance;
  double absolute_tolerance;
  std::size_t max_evaluations;
};

struct OpticalDepthResult {
  OpticalDepth optical_depth;
  std::size_t evaluations;
};

struct ConversionProbabilityResult {
  OpticalDepth optical_depth;
  Probability probability;
  std::size_t evaluations;
};

OpticalDepthResult optical_depth(
    Energy photon_energy,
    const StraightTrajectory& trajectory,
    const MagneticFieldModel& field,
    Length upper_path_length,
    const IntegrationOptions& options);

ConversionProbabilityResult conversion_probability(
    Energy photon_energy,
    const StraightTrajectory& trajectory,
    const MagneticFieldModel& field,
    const IntegrationOptions& options);
```

`conversion_probability()` calls `optical_depth()` at `trajectory.length` and
then computes `P = -expm1(-tau)`. Sampling calls `optical_depth()` first for
`tau_total`, then repeatedly for intermediate path lengths during bisection.

Sampling:

```cpp
struct LocalizationOptions {
  Length absolute_path_tolerance;
  double relative_path_tolerance;
  std::size_t max_iterations;
};

struct ConversionSample {
  bool converted;
  double random_uniform;
  OpticalDepth target_optical_depth;
  OpticalDepth total_optical_depth;
  Probability total_probability;
  std::optional<Length> path_length;
  std::optional<Position> position;
  std::size_t evaluations;
};

ConversionSample sample_first_conversion(
    Energy photon_energy,
    const StraightTrajectory& trajectory,
    const MagneticFieldModel& field,
    RandomEngine& rng,
    const IntegrationOptions& integration_options,
    const LocalizationOptions& localization_options);
```

`sample_first_conversion()` draws `U`, computes `tau_target = -log(U)`, compares
against `tau_total`, and if conversion occurs uses bisection to return the first
compatible upper bracket within the path tolerance.

## 4. Resultados de Referencia

For:

```text
E = 7e19 eV
B_perp = 2.115138e-5 T
L = 1e7 m
```

Current reference values:

```text
rate = 83.212547954233116 s^-1
tau = 2.7756718267486611
probability = 0.93769239728892728
```

Sampling references for the same geometry:

```text
seed 12345 converts at s ~= 3704534.2015335336 m
seed 43 does not convert
```

## 5. Comandos de Verificacion

Configure:

```sh
cmake -S . -B build
```

Build:

```sh
cmake --build build
```

Tests:

```sh
ctest --test-dir build --output-on-failure
```

Representative CLI commands:

```sh
./build/preshaires rate \
  --energy-eV 7e19 \
  --Bperp-T 2.115138e-5

./build/preshaires probability \
  --energy-eV 7e19 \
  --length-m 1e7 \
  --direction 0,0,-1 \
  --uniform-field-T 2.115138e-5,0,0

./build/preshaires sample-conversion \
  --energy-eV 7e19 \
  --length-m 1e7 \
  --direction 0,0,-1 \
  --uniform-field-T 2.115138e-5,0,0 \
  --seed 12345

./build/preshaires sample-conversion \
  --energy-eV 7e19 \
  --length-m 1e7 \
  --direction 0,0,-1 \
  --uniform-field-T 2.115138e-5,0,0 \
  --seed 43
```

## 6. Deuda Conocida

- Pair energy distribution.
- Creation of `e+` and `e-` particles.
- Segmentation by field-table nodes to better handle discontinuities or sharp
  changes.
- AIRES RNG adapter.
- Stable C API.
- Magnetic emission and cascade.

## 7. Siguiente Etapa del Roadmap

v0.4 should focus on the pair energy distribution.

Before implementing it, audit the physical formula used by MaGICS: primary
source, normalization, electron/positron symmetry and validity range. Do not
copy `pair_production_spectrum()` blindly. It is legacy C code coupled to
implicit conventions, global state and AIRES RNG usage.

## 8. Riesgos Abiertos

- The license of inherited MaGICS code remains unresolved.
- `src/NR` may carry Numerical Recipes restrictions.
- Do not reuse `src/NR` in the new core.
- Earth geometry and IGRF remain deferred.
