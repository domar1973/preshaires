# Erber Formula Comparison

This document compares the local MaGICS magnetic pair-conversion formula with the primary Erber 1966 convention available locally as `erber1966.pdf` (`High-Energy Electromagnetic Conversion Processes in Intense Magnetic Fields`, Rev. Mod. Phys. 38, 626, DOI in PDF metadata: `10.1103/RevModPhys.38.626`).

No production code is changed by this audit.

## A. Formula Implemented by MaGICS

From `src/MaGICS/geometry.c:195-197`, MaGICS computes:

```text
chi_M = (E_num * e / (2 m_e c^2)) * (Bperp_num / Bcrit)
```

where:

- `E_num` is numerically the particle energy in GeV.
- `Bperp_num` is numerically the transverse geomagnetic field in nT, because `geomagnetic_()` returns nT.
- The numerical product is nevertheless equal to using `E` in eV and `Bperp` in tesla:

```text
E_GeV * B_nT = (E_eV / 1e9) * (B_T * 1e9) = E_eV * B_T
```

Thus the physically transparent form is:

```text
chi_M = (1/2) * (E_eV * e / (m_e c^2)) * (Bperp_T / Bcrit)
```

From `src/MaGICS/functions.c:134-143`, MaGICS computes a linearized step probability:

```text
p_step = Gamma_M * dt
Gamma_M = 1.234e18 s^-1 * (m_e c^2 / E_J) * K_{1/3}(2 / (3 chi_M))^2
E_J = E_GeV * 1e9 * e
```

`Gamma_M` is a rate per unit time. The attenuation coefficient per unit length is:

```text
alpha_M = Gamma_M / c
```

Constants from `src/MaGICS/constants.h:3-8`:

| Constant | Code value | Unit |
| --- | ---: | --- |
| `c` | `299792458` | m/s |
| `hbar` | `1.05457e-34` | J s |
| `e` | `1.6021892e-19` | C, also used as J/eV |
| `electron_mass` | `9.109534e-31` | kg |
| `alpha` | `1/137.0` | dimensionless |

The critical field is computed at `src/MaGICS/init.c:169-170`:

```text
Bcrit = (m_e c)^2 / (e hbar)
```

The formula above is directly from code. Interpreting `Bcrit` as tesla and resolving the GeV/nT cancellation is interpretation supported by the diagnostic CSV and earlier local audit.

## B. Primary Erber Convention

The local primary source `erber1966.pdf` contains Sec. 3B, "Magnetic Pair Production." The page image for p. 642 shows:

- Eq. (3.1): pair creation over path length `d` is controlled by a photon attenuation coefficient `alpha(x)`, with `1 - exp[-alpha(x)d]` and the small-depth approximation `alpha(x)d`.
- The photon parameter is:

```text
x = (1/2) * (h nu / m c^2) * (H / Hcr)
```

- Eq. (3.3a):

```text
alpha(x) = (1/2) * (alpha / lambda_c) * (H / Hcr) * T(x)
```

- Eq. (3.3d), analytic approximation:

```text
T(x) ~= 0.16 * x^-1 * K_{1/3}(2 / (3 x))^2
```

- Eq. (3.4), after combining:

```text
alpha(x) = 0.16 * (alpha / lambda_c) * (m c^2 / h nu)
           * K_{1/3}( (4/3) * (m c^2 / h nu) * (Hcr / H) )^2
```

This is an attenuation coefficient per unit length. Multiplying by `c` gives a time rate:

```text
Gamma_E = c * alpha(x)
        = 0.16 * alpha * (m_e c^2 / hbar) * (m_e c^2 / E)
          * K_{1/3}(2 / (3 x))^2
```

This assumes the Erber setup on the cited page: photon propagation perpendicular to a uniform magnetic field, `H << Hcr`, ultrarelativistic limit, leading terms retained, magnetic and electron-positron bound-state resonances neglected. The page does not state a polarization split in this approximate formula; it is treated here as the total attenuation coefficient as presented by Erber.

## C. Convention Table

| Source / implementation | Quantum parameter | Definition | Factor `1/2` | Bessel argument | Prefactor coefficient | Dimensional scale | Time or length | Polarization | Approximation regime |
| --- | --- | --- | --- | --- | ---: | --- | --- | --- | --- |
| MaGICS | `chi_M` | `(1/2)(E/m_e c^2)(Bperp/Bcrit)` after GeV/nT cancellation | Present | `2/(3 chi_M)` | `C_M = 0.21775772588580714` | `alpha m_e c^2 / hbar` | time rate | not explicit | labeled Erber approximation in code |
| Erber 1966 Eq. (3.3a, 3.3d, 3.4) | `x` | `(1/2)(h nu/m c^2)(H/Hcr)` | Present | `2/(3 x)` | `C_E = 0.16` | length: `alpha/lambda_c`; time: `alpha m_e c^2/hbar` | length, convertible to time by `c` | not split on cited page | `H << Hcr`, ultrarelativistic, no bound-state resonances |
| MaGICS local comment | `chi` | comment says `(E/m)(B/Bcr)` | Absent in comment | not specified by comment | not specified | not specified | not specified | not specified | comment only; implementation differs |
| PRESHOWER | not used | no local formula found in this repo/AIRES tree | not available | not available | not available | not available | not available | not available | not used as authority |

## D. Algebraic Transformations

Let a source convention define:

```text
chi_E = a * chi_M
```

To keep the Bessel argument identical, a formula written as `K_{1/3}(2/(3 chi_E))^2` requires `a = 1`. If `a != 1`, the argument changes:

```text
K_{1/3}(2/(3 chi_E))^2 = K_{1/3}(2/(3 a chi_M))^2
```

No constant prefactor can make this equal to `K_{1/3}(2/(3 chi_M))^2` for all `chi_M`, because the Bessel function is nonlinear.

Specific cases:

### `chi_E = chi_M`

```text
K_{1/3}(2/(3 chi_E))^2 = K_{1/3}(2/(3 chi_M))^2
```

The Erber parameter `x` equals `chi_M` exactly when `E = h nu` and `Bperp = H`.

The full formulas then differ only by the dimensionless prefactor:

```text
Gamma_M / Gamma_E = C_M / C_E
```

### `chi_E = 2 chi_M`

```text
K_{1/3}(2/(3 chi_E))^2 = K_{1/3}(1/(3 chi_M))^2
```

This is not algebraically equivalent to MaGICS for all `chi_M`. Changing the prefactor cannot repair the changed Bessel argument globally.

### `chi_E = chi_M / 2`

```text
K_{1/3}(2/(3 chi_E))^2 = K_{1/3}(4/(3 chi_M))^2
```

Again, no constant prefactor restores equality for all `chi_M`.

### Erber Eq. (3.4) Argument

Erber's Eq. (3.4) argument:

```text
(4/3) * (m c^2 / E) * (Hcr / H)
```

Using `x = (1/2)(E/m c^2)(H/Hcr)`:

```text
2/(3x) = 2 / [3 * (1/2)(E/m c^2)(H/Hcr)]
       = (4/3)(m c^2/E)(Hcr/H)
```

So MaGICS' `2/(3 chi_M)` is exactly Erber's argument if `chi_M = x`.

## E. Prefactor Reconstruction

Using MaGICS constants:

```text
alpha * m_e c^2 / hbar = 5.666848305750188e18 s^-1
```

MaGICS:

```text
1.234e18 s^-1 = C_M * alpha * m_e c^2 / hbar
C_M = 0.21775772588580714
```

Erber Eq. (3.4), converted from length to time:

```text
Gamma_E = 0.16 * alpha * m_e c^2 / hbar
          * (m_e c^2 / E)
          * K_{1/3}(2/(3x))^2
```

Therefore:

```text
C_E = 0.16
C_M / C_E = 1.3609857867862947
```

The difference is not a notational change in `chi`: with `chi_M = x`, the Bessel argument and energy dependence are identical. It is a 36.1% prefactor difference relative to Erber's Eq. (3.4) approximation as read from the local primary source.

Relevant factors:

| Factor | MaGICS | Erber Eq. (3.4) time-rate form |
| --- | --- | --- |
| `alpha` | embedded in numeric `1.234e18` if interpreted as QED scale | explicit |
| `lambda_c` | not explicit | appears in length coefficient; `c/lambda_c = m c^2/hbar` |
| `hbar` | embedded in QED scale comparison and `Bcrit` | explicit via time conversion |
| `B/Bcrit` | inside `chi_M` | inside `x`; cancels out of Eq. (3.4)'s explicit prefactor via `T(x) ~ x^-1` |
| `E/(m c^2)` | inside `chi_M`; inverse appears in rate prefactor | same |
| numeric coefficient | `0.21775772588580714` | `0.16` |

## F. Numerical Comparison

The companion script `tools/compare_erber_conventions.py` evaluates:

- MaGICS with physical inputs `E_eV` and `Bperp_T`;
- Erber 1966 Eq. (3.4) converted to a time rate.

Because `chi_M = x` and the Bessel arguments are identical, the ratio is constant wherever both rates are finite:

```text
rate_MaGICS / rate_Erber1966 = 1.3609857867862947
```

The diagnostic case requested:

```text
E = 7e19 eV
Bperp = 2.115138e-5 T
```

is included in `diagnostics/erber_convention_comparison.csv`.

## G. Classification

### Demonstrated

- MaGICS implements `chi_M = (1/2)(E/m_e c^2)(Bperp/Bcrit)` numerically, once the GeV/nT cancellation is made explicit.
- Erber 1966 defines `x = (1/2)(h nu/m c^2)(H/Hcr)`.
- Therefore `chi_M = x` for a photon of energy `E = h nu` and transverse field `Bperp = H`.
- MaGICS and Erber use the same Bessel argument under that identification.
- MaGICS' prefactor corresponds to `C_M = 0.21775772588580714` multiplying `alpha m_e c^2/hbar`.
- Erber Eq. (3.4), converted to time rate, uses `C_E = 0.16`.
- The prefactor ratio is `1.3609857867862947`.

### Inferido

- MaGICS likely intended to implement Erber's analytic approximation because the code comment says "Erber's analytical aproximation" and the argument matches Erber's Eq. (3.3d)/(3.4).
- The use of `Bperp` generalizes Erber's perpendicular-propagation setup by replacing `H` with the field component transverse to photon motion.

### No Resuelto

- The local MaGICS repository does not explain the origin of the numeric coefficient `1.234e18`.
- It is not resolved whether MaGICS intentionally used a different polarization average, correction factor, or another source convention that changes `0.16` to `0.21775772588580714`.
- The exact relation to Homola/Rygielski or PRESHOWER was not established locally because no accessible local formula was found in this tree.

### Posibles Defectos

- Relative to Erber 1966 Eq. (3.4) as the primary convention, MaGICS is not exactly equivalent: it is larger by the constant prefactor ratio `1.3609857867862947`.
- This is a demonstrated algebraic/numerical discrepancy with that equation, not a proven bug in MaGICS, because the source of MaGICS' coefficient remains undocumented.
- The GeV/nT cancellation in `chi()` remains a fragile implementation convention, but it does not explain the Erber prefactor discrepancy.
