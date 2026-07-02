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

The first bootstrap formula is deterministic and uses no RNG. Future stochastic
transport must receive RNG state through explicit parameters or interfaces.
There must be no mutable global RNG in the core.

## I/O

The core performs no file or console I/O. The CLI prints machine-readable text
for the current diagnostic command. Future I/O formats should be introduced only
after the core data model is stable.

## Tests

The first unit tests cover unit conversions, `chi`, the Erber rate, finiteness,
positivity and local monotonicity. Historical diagnostics under `docs/`,
`diagnostics/` and `tests/` remain reference material for later migration.

## Migration Strategy

1. Bootstrap the C++20 core with the Erber local rate and explicit units.
2. Add focused reference tests around each migrated formula.
3. Introduce transport by optical-depth sampling rather than fixed linearized
   step probabilities.
4. Move AIRES coupling into `adapters/aires` without leaking AIRES into core.
5. Preserve legacy MaGICS until equivalent behavior and tests exist.
