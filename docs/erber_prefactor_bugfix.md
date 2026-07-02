# Erber Pair-Production Prefactor Bugfix

## Summary

`pair_production_probability()` used a hard-coded dimensional prefactor:

```c
1.234E18 * deltatime * (electron_mass * c * c / (Energy * 1E9 * e))
    * K_{1/3}(2 / (3 chi))^2
```

The corrected expression is:

```c
0.16 * alpha * electron_mass * c * c / hbar
    * deltatime * (electron_mass * c * c / (Energy * 1E9 * e))
    * K_{1/3}(2 / (3 chi))^2
```

The change replaces an opaque precomputed rate scale with the Erber 1966
coefficient and the existing MaGICS physical constants.

## Dimensional Derivation

Erber Eq. (3.4), using the approximation in Eq. (3.3d), gives the attenuation
coefficient by length. Converting it to a rate by time multiplies by `c`:

```text
Gamma =
    0.16
    * alpha
    * m_e c^2 / hbar
    * m_e c^2 / E
    * K_{1/3}(2 / (3 chi))^2
```

The factor `0.16` is dimensionless. `alpha` and the two energy ratios are
dimensionless except for `m_e c^2 / hbar`, which has units `s^-1`. Therefore
`Gamma` is a rate per unit time, and the per-step probability remains
`Gamma * deltatime`.

Using the MaGICS constants:

```text
0.16 * alpha * m_e c^2 / hbar = 9.066957289200302e17 s^-1
```

## Previous Ratio

The old literal was:

```text
1.234e18 s^-1
```

The exact ratio with MaGICS constants is:

```text
1.234e18 / (0.16 * alpha * m_e c^2 / hbar)
    = 1.3609857867862947
```

Thus the old local rate was larger than the Erber rate by about `36.1%` for
any row where the Bessel value is finite and all other inputs are unchanged.

## Impact

The local pair-production rate is reduced by the constant factor:

```text
1 / 1.3609857867862947 = 0.7347615242463777
```

The integrated conversion probability is not reduced by that exact factor in
general. MaGICS accumulates survival over many steps, so the final probability
depends nonlinearly on the step probabilities:

```text
P = 1 - product_i(1 - p_i)
```

For small total optical depth this is close to a linear scaling. At larger
depths, saturation makes the final probability respond less than linearly.

## Scope

Changed:

- `pair_production_probability()` now builds the rate prefactor from
  `0.16 * alpha * electron_mass * c * c / hbar`.
- The independent Python audit uses the same corrected Erber prefactor.
- Regression tests cover the prefactor and the documented `E = 7e19 eV`,
  `B_perp = 2.115138e-5 T`, `chi = 0.3281999188925183` case.

Explicitly not changed:

- the current `chi` definition;
- the existing numerical GeV/nT cancellation in `chi()`;
- the Bessel argument `2 / (3 chi)`;
- `dbskr3_()`;
- the integration scheme;
- random generation;
- `inic()` dispatch;
- AIRES;
- existing audit documentation.

## Historical Note

The previous `1.234E18` literal was a manual prefactor calculation error made
approximately 23 years ago, as confirmed by the original MaGICS author. The
bug was algebraic only in the prefactor: the demonstrated MaGICS `chi`
definition equals Erber's photon parameter `x`, and the Bessel argument and
energy-field dependence match Erber's convention.

## Primary Reference

Thomas Erber, "High-Energy Electromagnetic Conversion Processes in Intense
Magnetic Fields", Reviews of Modern Physics 38, 626 (1966), especially
Eqs. (3.3d) and (3.4).
