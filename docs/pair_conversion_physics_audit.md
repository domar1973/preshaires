# Pair Conversion Physics Audit

This note audits the local conversion-rate formula as a coupled set:

```text
chi = (E_GeV * e / (2 m_e c^2)) * (Bperp_nT / Bcrit)
rate = 1.234e18 * (m_e c^2 / E_J) * K_{1/3}(2 / (3 chi))^2
E_J = E_GeV * 1e9 * e
```

No production formula is changed here.

## Demonstrated Locally

- `pair_production_probability()` implements `1.234E18 * dt * (m_e c^2 / (Energy * 1E9 * e)) * K_{1/3}(2/(3 chi))^2` in `src/MaGICS/functions.c:135-143`.
- `chi()` implements `(this_particle->energy * e / (2 * electron_mass * c * c)) * (bperp / bcrit)` in `src/MaGICS/geometry.c:195-197`.
- `bcrit` is computed from SI constants as `(electron_mass * c)^2 / (e * hbar)` in `src/MaGICS/init.c:169-170`.
- `geomagnetic_()` returns magnetic-field strength numerically in nT in the diagnostic run: AIRES reports `22.665 uT`; the CSV reports `B = 22661.8`; `22661.8 nT = 22.6618 uT`.
- `dbskr3_(x, 1)` returns the unscaled modified Bessel function `K_{1/3}(x)`: `src/NR/dbskr3_fromNR.f:39-42` sets `nu = abs(nu3)/3`, calls `bessik`, and returns `rk`; `src/NR/bessik.f:25-30` documents `rk = K_nu`.

## Why the `1e9` Cancellation Works

The physically transparent SI expression would use:

```text
E_J = E_GeV * 1e9 * e
B_T = B_nT * 1e-9
```

The current `chi()` instead uses the numeric value `E_GeV` where an eV-like value would be expected by the multiplication `Energy * e`, and it also uses the numeric value `B_nT` where a tesla value would be expected by division by SI-valued `Bcrit`.

The product contains both factors:

```text
E_GeV * B_nT = (E_eV / 1e9) * (B_T * 1e9) = E_eV * B_T
```

Thus the two missing conversions cancel exactly in the product that enters `chi`.

## Why This Is Fragile

The cancellation is dimensional bookkeeping hidden in two different subsystems:

- particle energy is carried numerically as GeV;
- AIRES `geomagnetic_()` returns the magnetic field numerically as nT;
- `bcrit` is computed from SI constants;
- the formula in `chi()` contains no explicit `1e9` or `1e-9`.

If someone corrected only the energy conversion in `chi()` by replacing `Energy * e` with `Energy * 1e9 * e`, `chi` would become `1e9` times larger. If someone corrected only the magnetic field by converting `Bperp_nT` to tesla, `chi` would become `1e9` times smaller. Either isolated change would severely alter the rate because the Bessel argument is `2/(3 chi)`.

Any unit cleanup must convert energy and magnetic field together, with tests that preserve the current numeric `chi` for known diagnostic rows.

## Photon Quantum Parameter Convention

The local comment in `src/MaGICS/geometry.c:136-138` says `chi=(E/m)*(B/B_cr)`, but the implemented formula includes a factor `1/2`:

```text
chi_code = (E / (2 m_e c^2)) * (Bperp / Bcrit)
```

The pair-production rate then uses `K_{1/3}(2/(3 chi_code))^2`.

The repository identifies this as "Erber's analytical approximation" in `src/MaGICS/functions.c:134`, but it does not include the primary Erber equation. Therefore the following is unresolved from local evidence alone:

- whether the code's `chi` is the standard photon quantum parameter or one half of it;
- whether `2/(3 chi)` is paired with the code's `chi` definition or with a differently normalized parameter in the source formula;
- whether the coefficient `1.234e18` was chosen for this exact `chi` convention.

The primary reference needed is the original Erber pair-production attenuation-rate equation, including its definition of the photon field parameter and Bessel argument.

## Prefactor Decomposition

Using the constants in `src/MaGICS/constants.h:3-8`:

```text
alpha * m_e c^2 / hbar = 5.666848305750188e18 s^-1
```

The coded prefactor can be written as:

```text
1.234e18 s^-1 = 0.21775772588580714 * alpha * m_e c^2 / hbar
```

This demonstrates the natural QED rate scale and the dimensionless coefficient implied by MaGICS' constants. The repository does not provide a derivation for the coefficient `0.21775772588580714`.

## Compatibility Assessment

Numerically, the CSV audit separates two facts:

1. The implemented C/Fortran rate is internally reproducible by an independent Python evaluation using `scipy.special.kv(1/3, x)`.
2. The dimensional convention in `chi()` is implicit and fragile but currently cancels to the same numeric value as `E_eV * B_T`.

Physically, local sources are insufficient to prove that the prefactor, the `1/2` in `chi`, and the argument `2/(3 chi)` all belong to the same published convention. The code may be self-consistent numerically while still requiring a primary-reference check for the normalization of `chi` and the attenuation-rate prefactor.

## Not Resolved Here

- No correction is made to `pair_production_probability()`.
- No correction is made to `chi()`.
- No conversion of AIRES nT to tesla is introduced.
- No dispatch issue in `inic()` is corrected.
