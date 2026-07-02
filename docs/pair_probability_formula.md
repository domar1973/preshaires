# `pair_production_probability()` numerical audit

## Implemented formula

The implementation is in `src/MaGICS/functions.c:127-136`:

```c
besselarg=2/(3*chi);
result=1.234E18*deltatime*(electron_mass*c*c/(Energy*1E9*e))*
             sqr(dbskr3_(&besselarg,&mu));
```

With `mu = 1`, `dbskr3_(&x,&mu)` returns `K_{1/3}(x)`. The exact formula used by the code is therefore:

```text
P_step = 1.234e18 * dt * (m_e c^2 / E_gamma) * K_{1/3}(2 / (3 chi))^2
```

where `E_gamma = Energy * 1e9 * e` joules because `Energy` is passed in GeV.

## Units

| Factor | Code symbol | Unit in MaGICS | Evidence |
| --- | --- | --- | --- |
| `1.234e18` | literal | `s^-1` | Must cancel `dt`; no declaration elsewhere. |
| `dt` | `deltatime` | `s` | `conversion_probability()` propagates with `c * dt`; `c` is m/s in `src/MaGICS/constants.h:3`. |
| `m_e` | `electron_mass` | `kg` | `src/MaGICS/constants.h:7`. |
| `c` | `c` | `m s^-1` | `src/MaGICS/constants.h:3`. |
| `m_e c^2` | `electron_mass*c*c` | `J` | SI constants. |
| `Energy` | argument | `GeV` | Converted by `Energy * 1E9 * e` in `src/MaGICS/functions.c:134`. |
| `e` | `e` | `C`, used here as `J/eV` | `src/MaGICS/constants.h:5`; the numeric value is the elementary charge. |
| `chi` | argument | dimensionless | Computed as an energy ratio times a field ratio in `src/MaGICS/geometry.c:196`; see the implicit unit convention below. |
| `K_{1/3}` | `dbskr3_` | dimensionless | Function of dimensionless `2/(3 chi)`. |

The returned value is dimensionless.

## Prefactor `1.234E18`

`1.234E18` appears only as a hard-coded literal in `src/MaGICS/functions.c:134`. In dimensional terms it is a rate prefactor in `s^-1`.

The nearest standard scale available from the constants declared by MaGICS is:

```text
alpha * m_e c^2 / hbar = 5.666848305750188e18 s^-1
```

using `alpha`, `electron_mass`, `c`, and `hbar` from `src/MaGICS/constants.h:3-8`. The coded value is:

```text
1.234e18 / (alpha * m_e c^2 / hbar) = 0.21775772588580714
```

No expression, comment, or named constant in this repository derives `1.234E18`; its direct origin in MaGICS is the literal in `pair_production_probability()`. The surrounding comment labels the approximation as Erber's analytical approximation at `src/MaGICS/functions.c:127`.

## Exact `chi` used by the code

`conversion_probability()` calls `chi(auxiliar_photon)` at `src/MaGICS/functions.c:225`.

`chi()` is implemented in `src/MaGICS/geometry.c:139-197`. It:

1. Rotates particle position and direction from AIRES coordinates to geocentric Cartesian coordinates.
2. Computes geographic latitude, longitude and altitude.
3. Calls `cartesian_magnetic_field()`, which calls AIRES `geomagnetic_()` and returns the geomagnetic field numerically in nanotesla.
4. Computes the field transverse to the photon direction:

```text
Bperp = |B - k (B . k)|
```

5. Returns:

```text
chi = (Energy * e / (2 * electron_mass * c^2)) * (Bperp / bcrit)
```

This is exactly `src/MaGICS/geometry.c:196`. The convention is implicit:

- `particle.energy` is numerically in GeV.
- `geomagnetic_()` returns `B` numerically in nT.
- The code multiplies `Energy` by `e`, not by `1e9 * e`, and divides the nT-valued `Bperp` by the SI-valued numeric `bcrit`.
- The missing `1e9` in the energy conversion and the `1e-9` from using nT instead of T cancel numerically.

Therefore the current `chi()` value is numerically correct, but the unit convention is fragile. One must not "fix" only the energy unit or only the magnetic-field unit: changing just one side would change `chi` by a factor of `1e9`.

`bcrit` is set in `inic()` as:

```text
bcrit = (electron_mass * c)^2 / (e * hbar)
```

at `src/MaGICS/init.c:169-170`, in tesla under the code's SI convention.

## Rate or probability?

`pair_production_probability()` returns a linearly integrated per-step probability:

```text
P_step = rate(chi, E) * dt
```

It is not an exact finite-step probability of the form `1 - exp(-rate * dt)`.

`conversion_probability()` then multiplies survival factors:

```text
survival *= (1 - P_step)
return 1 - survival
```

at `src/MaGICS/functions.c:228-235`. Therefore the integrated result is exact only for the discrete product of the linearized step probabilities; it is not the continuous exponential survival integral unless all step probabilities are small.

## When can a step probability exceed 1?

The code has no clamp. A step can exceed 1 whenever:

```text
1.234e18 * dt * (m_e c^2 / E_gamma) * K_{1/3}(2/(3 chi))^2 > 1
```

This can happen for sufficiently large `dt`, sufficiently low `Energy`, or sufficiently large `chi` because the Bessel factor grows as the argument `2/(3 chi)` becomes small. `conversion_probability()` uses `dt = 1e-5 s` above `0.1 * earth_radius` and `dt = 5e-7 s` below that altitude at `src/MaGICS/functions.c:225-228`, so oversized step probabilities are possible if the local rate exceeds `1/dt`.

If `P_step > 1`, the survival update `survival *= (1 - P_step)` changes sign and the accumulated probability can become unphysical. If `P_step < 0`, which the formula should not produce for positive finite inputs but can occur through invalid/NaN-adjacent numerical states, the survival factor exceeds one.
