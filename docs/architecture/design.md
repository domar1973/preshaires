# preshaires Architecture Bootstrap

## Objectives

- Provide a C++20 physical core shared by standalone and AIRES-coupled modes.
- Make units explicit at public boundaries.
- Keep formulas auditable through named constants, references, equations,
  assumptions and validity ranges.
- Preserve historical MaGICS behavior as reference material while migrating
  progressively.
- Favor correctness, clarity, reproducibility and tests over premature
  optimization.

## Non-Objectives

- No full transport rewrite in this stage.
- No real AIRES linkage in this stage.
- No Python bindings, YAML configuration, plugin mechanism or field-model zoo.
- No deletion or relocation of MaGICS legacy code.

## Layers

1. `preshaires-core`: C++20 physical formulas and future transport primitives.
2. Stable C API: planned, not implemented in this bootstrap.
3. `preshaires-cli`: standalone command-line access to core functionality.
4. `preshaires-aires`: adapter layer under `adapters/aires` only.
5. Python reference: existing scripts remain reference material; bindings are
   deferred.

## Allowed Dependencies

- `preshaires-core` may depend on the C++20 standard library only at this stage.
- `preshaires-cli` may depend on `preshaires-core`.
- `adapters/aires` may depend on `preshaires-core` and, later, AIRES.
- AIRES headers and symbols are forbidden outside `adapters/aires`.

## Standalone / AIRES Separation

The standalone CLI owns parsing and presentation only. AIRES integration will
own conversion between AIRES event data and typed core inputs. Neither layer
owns physical formulas. The core receives explicit values such as photon energy
in eV and transverse magnetic field in tesla.

## Units Model

The bootstrap uses small explicit wrappers:

- `Energy { double eV; }`
- `MagneticField { double tesla; }`
- `Length { double meter; }`
- `Time { double second; }`
- `Rate { double per_second; }`

This avoids a dependency before the required API shape is clear. If the type
surface grows enough to justify a units library, that decision should be made
with migration cost, compile time and API stability documented.

The named physical constants currently use the same numerical values as the
MaGICS diagnostic convention so the historical reference case remains exactly
reproducible. A later constants update must be done deliberately and with
reference-test changes.

## Error Model

The current deterministic formula returns `NaN` for invalid nonpositive energy
or nonfinite input, and zero rate for zero transverse field. A richer error
model should be introduced before broadening the public C API.

## RNG

Stochastic operations receive RNG state through `RandomEngine`, whose
`uniform_open01()` contract is a finite value `0 < U < 1`. There is no mutable
global RNG in the core. Standalone CLI code may adapt `std::mt19937_64`; AIRES
can later provide its own adapter through the same interface.

## I/O

The core performs no file or console I/O. The CLI prints machine-readable text
for the current diagnostic command. Future I/O formats should be introduced only
after the core data model is stable.

For v0.2, tabulated magnetic-field input is read only by the CLI. The core
receives already parsed C++ nodes with explicit units.

## Straight-Trajectory Model

v0.2 adds a prescribed rectilinear trajectory:

- `Position` stores Cartesian coordinates in meters.
- `Direction` is a unit vector built with `make_direction()`.
- `StraightTrajectory` stores a start position, direction and finite
  nonnegative path length.

This is not Earth geometry. It is a minimal integration domain for optical
depth tests and standalone diagnostics.

## Magnetic-Field Interface

`MagneticFieldModel` supplies `field_at(position, path_length)`. Two concrete
models are available:

- `UniformMagneticField`, which returns one fixed vector.
- `TabulatedMagneticField`, which linearly interpolates strictly increasing
  path-length nodes and refuses extrapolation.

The transverse field is computed as `|B - (B dot n)n|`, with small negative
roundoff residues clamped to zero.

## Conversion-Probability Integration

The deterministic probability is:

```text
tau = integral_0^L [Gamma(E, B_perp(s)) / c] ds
P = -expm1(-tau)
```

`Gamma` is provided by the existing Erber local-rate function; the integration
code does not duplicate the pair-production formula. The numerical method is an
iterative adaptive Simpson rule with explicit absolute/relative tolerances,
evaluation counting and a maximum-evaluation failure path. It computes only
survival probability, not the sampled interaction point.

## First-Conversion Sampling

v0.3 adds optical-depth sampling:

```text
U in (0, 1)
tau_target = -log(U)
```

If `tau_target < tau_total`, the first conversion point is localized by
bisection on `tau(s) - tau_target` using repeated calls to the same
`optical_depth()` implementation. The returned point is the first compatible
upper bracket within the requested path-length tolerance. This is not fixed-step
Monte Carlo transport and does not create pair particles.

## Deferred Geometry and IGRF

Earth coordinates, altitude, AIRES coordinate transforms and IGRF remain
deferred because they would couple this stage to site/date conventions and
field-model ownership. v0.2 instead establishes the field-provider interface
that those later systems can implement.

## Tests

The first unit tests cover unit conversions, `chi`, the Erber rate, finiteness,
positivity and local monotonicity. Historical diagnostics under `docs/`,
`diagnostics/` and `tests/` remain reference material for later migration.

v0.2 extends tests to geometry, uniform and tabulated fields, analytic uniform
optical depth, parallel-field zero probability, linear-profile integration,
convergence, small probabilities and saturated probabilities.

v0.3 adds deterministic RNG tests, invalid-RNG rejection, analytic uniform
sampling, no-conversion cases, zero-rate/plateau localization, reproducibility
and a seeded statistical conversion-fraction test.

## Migration Strategy

1. Bootstrap the C++20 core with the Erber local rate and explicit units.
2. Add focused reference tests around each migrated formula.
3. Introduce transport by optical-depth sampling rather than fixed linearized
   step probabilities.
4. Move AIRES coupling into `adapters/aires` without leaking AIRES into core.
5. Preserve legacy MaGICS until equivalent behavior and tests exist.
