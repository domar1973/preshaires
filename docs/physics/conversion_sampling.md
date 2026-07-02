# Sampling the First Conversion

## Scope

v0.3 samples only whether and where the first geomagnetic pair conversion
occurs along an already prescribed straight trajectory. It does not create the
electron-positron pair, sample the pair energy split, emit magnetic radiation,
run a cascade, evaluate Earth geometry, call IGRF or couple to AIRES.

## Optical-Depth Target

For one photon, draw:

```text
U ~ Uniform(0, 1)
tau_target = -log(U)
```

The RNG contract is strict: the core accepts only finite values with
`0 < U < 1`. Invalid values are errors, not clipped values.

## Conversion Criterion

Let:

```text
tau_total = integral_0^L [Gamma(E, B_perp(s)) / c] ds
```

Then:

- if `tau_total = 0`, the photon never converts;
- if `tau_target >= tau_total`, the photon reaches the end of the trajectory;
- if `tau_target < tau_total`, the first conversion occurs at the first path
  length `s*` where `tau(s*) >= tau_target`.

The total conversion probability remains:

```text
P = -expm1(-tau_total)
```

## First-Point Distribution

For a constant attenuation coefficient
`alpha = Gamma(E, B_perp) / c`:

```text
tau(s) = alpha s
s* = -log(U) / alpha
```

The conversion before a finite endpoint `L` is accepted only when
`s* < L`, equivalent to `tau_target < alpha L`.

## Numerical Localization

The implementation first computes `tau_total` with the v0.2 adaptive Simpson
integrator. If conversion occurs, it brackets the root of:

```text
F(s) = tau(s) - tau_target
```

on `[0, L]` and uses bisection. Since `tau(s)` is monotone nondecreasing, the
upper bracket is always a compatible point. The returned path length is the
upper bracket after the interval width satisfies:

```text
high - low <= absolute_path_tolerance
              + relative_path_tolerance * max(|low|, |high|)
```

Using the upper bracket preserves the "first compatible point" convention up
to the requested path tolerance.

## Zero-Rate Intervals and Plateaus

Intervals with `B_perp = 0` have zero local rate and therefore flat optical
depth. If the target optical depth is reached at the beginning of a plateau,
bisection returns the beginning of that plateau within the path tolerance, not
the later end of the flat interval.

## Difference from Linear Stepping

The sampler does not use fixed Monte Carlo steps and does not approximate an
interaction probability as `Gamma * dt`. It samples a target optical depth and
then solves for the first point where the integrated optical depth reaches that
target.

## Limits of v0.3

- No pair creation.
- No pair energy distribution.
- No magnetic emission.
- No cascade.
- No Earth geometry or IGRF.
- No AIRES RNG adapter yet, although AIRES can later provide a `RandomEngine`.
