# Conversion Probability Along a Prescribed Path

## Scope

This note describes v0.2 only: deterministic conversion probability for a
photon of fixed energy moving along a prescribed straight trajectory. It does
not sample the conversion point, create a pair, propagate secondaries, evaluate
Earth geometry or call IGRF.

## Optical Depth

The local Erber pair-production function returns a time rate:

```text
Gamma(E, B_perp)  [s^-1]
```

For a path-length coordinate `s` in meters, the optical-depth increment is:

```text
d tau = [Gamma(E, B_perp(s)) / c] ds
```

Therefore:

```text
tau = integral_0^L [Gamma(E, B_perp(s)) / c] ds
```

The photon energy is fixed over this integral in v0.2.

## Units

- `Energy` is stored in eV.
- `Length` and path coordinate `s` are stored in meters.
- `MagneticFieldVector` components are stored in tesla.
- `Rate` is stored in `s^-1`.
- `Gamma / c` has units `m^-1`.
- `tau` and `P` are dimensionless.

## Transverse Field

For magnetic-field vector `B` and unit direction `n`:

```text
B_perp = |B - (B dot n) n|
```

The implementation evaluates the equivalent norm:

```text
B_perp^2 = |B|^2 - (B dot n)^2
```

Small negative values from floating-point roundoff are clamped to zero. Larger
negative values are treated as invalid numerical states.

## Survival Probability

The survival probability is `exp(-tau)`, so the conversion probability is:

```text
P = 1 - exp(-tau)
```

The code evaluates this as:

```text
P = -expm1(-tau)
```

This preserves precision when `tau` is small.

## Numerical Method

The integral is evaluated with an iterative adaptive Simpson rule. The method:

- uses the existing Erber local-rate function;
- applies absolute and relative tolerances;
- counts integrand evaluations;
- stops with an explicit error if `max_evaluations` is exceeded;
- uses an explicit stack instead of unbounded recursion;
- rejects invalid energy, trajectory length, tolerances and nonfinite
  integrand values.

## Analytic Test Cases

Uniform perpendicular field:

```text
tau = Gamma(E, B_perp) L / c
P = -expm1(-tau)
```

Uniform field parallel to the direction:

```text
B_perp = 0
tau = 0
P = 0
```

Linearly tabulated field:

```text
B(s) = B0 + (B1 - B0) s / L
```

This has no closed form for the Erber rate in the tests, so the adaptive result
is compared against a high-resolution independent trapezoidal reference built
inside the test executable.

## Limits of v0.2

- No Monte Carlo sampling.
- No interaction point.
- No electron-positron energy split.
- No magnetic emission.
- No atmospheric cascade.
- No Earth curvature or altitude model.
- No IGRF.
- No AIRES coupling.
