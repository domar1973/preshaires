# MaGICS magnetic conversion gamma -> e+ e- call graph

This document is a static reading of the repository as of the current
workspace state. It describes the code path that computes and samples magnetic
pair conversion of a gamma ray in MaGICS_TdF. No functional code changes are
implied here.

## 1. Entry point from AIRES to MaGICS

AIRES runs MaGICS as a special primary module. The example input registers the
binary with `SetGlobal Dyn MaGICS ./magics` and then defines special primaries
with `AddSpecialParticle ... {MaGICS}` in `example.inp:16`,
`example.inp:21`, `example.inp:25`, and `example.inp:32`.

The compiled executable is `magics`, linked against AIRES with `-lAires` in
`MakeMaGICS:162-165`. The C entry point is `main()` in
`src/MaGICS/magics.c:39`.

The AIRES handoff starts at `speistart_()`:

```text
AIRES special primary
  -> magics main()
     -> speistart_(&shower_number, &primary_energy,
                   default_injection_position, injection_depth,
                   ground_altitude, ground_depth,
                   d_ground_inj, shower_axis)
     -> inic()
     -> add_particle(PrimaryParticle, primary_energy, initial_position, ...)
     -> conversion_probability()
     -> propagate_particle()
     -> spaddp0_() for particles reaching the atmosphere
     -> speiend_()
```

Relevant lines:

- AIRES startup data are read in `src/MaGICS/magics.c:64-67`.
- One-time initialization is called in `src/MaGICS/magics.c:71-82`.
- The initial primary is inserted in `src/MaGICS/magics.c:99-103`.
- Conversion probability is computed and returned to AIRES variable 1 in
  `src/MaGICS/magics.c:105-112`.
- The Monte Carlo propagation loop calls the selected `propagate_particle`
  function in `src/MaGICS/magics.c:120-123`.
- Particles reaching atmosphere are sent back to AIRES by `spaddp0_()` in
  `src/MaGICS/magics.c:136-151`.
- Number of forced trials is returned to AIRES variable 2 in
  `src/MaGICS/magics.c:168-174`.

## 2. Full conversion probability call graph

The deterministic probability reported to AIRES is computed before the Monte
Carlo transport:

```text
main()
  -> inic()
     -> speigetparsc()
     -> croinputdata0_()
     -> getglobalc("PEMF_Deblevel")
     -> getglobalc("PEMF_init_altit")
     -> set bcrit
     -> set initial_position and propagate_particle
  -> add_particle(gamma)
  -> conversion_probability()
     -> clone first_particle into auxiliar_photon
     -> loop while altitude >= atm_injection_altitude
        -> chi(auxiliar_photon)
           -> rotAIRES2Cartesian(position)
           -> spheric_coords()
           -> heighfactor()
           -> cartesian_magnetic_field()
              -> geomagnetic_()
           -> rotCartesian2AIRES(B)
           -> compute B_perp
           -> return chi
        -> pair_production_probability(Energy, chi, dt)
           -> dbskr3_(2/(3 chi), 1)
        -> prob_non_up_to_now *= (1 - proba)
        -> position += c dt shower_axis
     -> return 1 - prob_non_up_to_now
  -> speisetrealvar_(1, probability)
```

Line anchors:

- `conversion_probability()` is in `src/MaGICS/functions.c:203-236`.
- `chi()` is in `src/MaGICS/geometry.c:139-198`.
- `cartesian_magnetic_field()` is in `src/MaGICS/geometry.c:60-87`.
- `pair_production_probability()` is in `src/MaGICS/functions.c:128-137`.
- `dbskr3_()` is declared in `src/MaGICS/prototypes.h:25` and implemented in
  `src/NR/dbskr3_fromNR.f`.

There is also `french_conversion_probability()` in
`src/MaGICS/functions.c:238-271`, using `frenchi()` in
`src/MaGICS/geometry.c:206-266`, but this path is not called from `main()`.
Despite the comment saying dipolar approximation, `frenchi()` currently calls
`cartesian_magnetic_field()` rather than `geomagneticr_()` at
`src/MaGICS/geometry.c:236-239`.

## 3. Functions by responsibility

| Responsibility | Function and lines | Notes |
| --- | --- | --- |
| Initial trajectory setup | `inic()` in `src/MaGICS/init.c:172-200` | Computes the initial position at altitude `PEMF_init_altit` by solving for `lambda` along `shower_axis`. |
| Time evolution / trajectory | `conversion_probability()` in `src/MaGICS/functions.c:221-231` | Straight-line propagation: `position[i] += c * dt * first_particle->motion_direction[i]`. Uses the initial gamma direction, not a curved path. |
| Monte Carlo gamma trajectory | `propagate_1D_klepikov_v2()` in `src/MaGICS/1D_klepikov_v2.c:21-77` and analogous functions | Straight-line steps until conversion or atmosphere. |
| Geomagnetic field | `cartesian_magnetic_field()` in `src/MaGICS/geometry.c:60-87` | Calls AIRES `geomagnetic_(&theta, &phi, &altitude, &fdate, &b, &inc, &dec)`. |
| Coordinate conversion | `rotAIRES2Cartesian()` in `src/MaGICS/geometry.c:94-113`, `rotCartesian2AIRES()` in `src/MaGICS/geometry.c:115-133`, `spheric_coords()` in `src/MaGICS/geometry.c:34-53` | Converts AIRES local coordinates to geocentric Cartesian and geographic coordinates. |
| Altitude | `heighfactor()` in `src/MaGICS/geometry.c:2-18`, used in `chi()` at `src/MaGICS/geometry.c:166` | Uses a series for `sqrt(1+x)-1` when `x < 0.6`. |
| `B_perp` | `chi()` in `src/MaGICS/geometry.c:189-193` | Projects out the component parallel to the photon direction: `B_perp = |B - k (B.k)|`. |
| `chi` | `chi()` in `src/MaGICS/geometry.c:195-197` | `chi = (E e / (2 m_e c^2)) * (B_perp / Bcrit)`. |
| Local conversion probability | `pair_production_probability()` in `src/MaGICS/functions.c:128-137` | Erber analytical approximation using `K_{1/3}` through `dbskr3_()`. |
| Integrated conversion probability | `conversion_probability()` in `src/MaGICS/functions.c:221-235` | Multiplies survival factors `(1 - proba)` over fixed time steps. |
| Conversion sampling | `propagate_1D_klepikov_v2()` in `src/MaGICS/1D_klepikov_v2.c:25-67`, `propagate_1D_klepikov()` in `src/MaGICS/1D_klepikov.c:13-48`, `propagate_1D_erberaprox_v2()` in `src/MaGICS/1D_erber_v2.c:25-59`, `propagate_1D_erberaprox()` in `src/MaGICS/1D_erber.c:15-48`, `propagate_3D_klepikov_v2()` in `src/MaGICS/others.c:49-91` | Draws `urandom_() < proba`; creates electron and positron with sampled or equal energy split. |
| Electron energy fraction | `random_electron_energy_fraction()` in `src/MaGICS/functions.c:101-125` | Builds a 1025-bin cumulative distribution from `pair_production_spectrum()`. |
| Pair spectrum | `pair_production_spectrum()` in `src/MaGICS/functions.c:88-97` | Klepikov spectrum: `K_{2/3}(alpha)/(v(1-v)) - kappa(alpha)/alpha`. |

## 4. Units used

The code mixes SI quantities with AIRES energies in GeV:

- `primary_energy`, `particle.energy`, and function argument `Energy` are in
  GeV. Evidence: `pair_production_probability()` converts `Energy` with
  `Energy * 1E9 * e` in `src/MaGICS/functions.c:134`.
- Positions, altitudes, distances, `earth_radius`, and `atmospheric_thickness`
  are in meters. Constants are declared in `src/MaGICS/constants.h:3-10`, and
  `read.me:117-120` states `PEMF_init_altit` is in meters.
- Time variables `mytime`, `last_time`, `dt`, `deltat`, and `emission_time`
  are in seconds. Propagation uses `c * dt` with `c` in m/s.
- Magnetic field `b`, `bx`, `by`, `bz`, and `bperp` are in tesla if AIRES
  `geomagnetic_()` returns SI field strength. This is required by the
  dimension of `bcrit`, computed in SI at `src/MaGICS/init.c:169-170`.
- `lastbfield.strength` stores `B/Bcrit`, not tesla, after
  `src/MaGICS/geometry.c:179-184`.
- `lastbfield.direction` is a dimensionless unit vector in AIRES coordinates,
  normalized in `src/MaGICS/geometry.c:182-183`.
- Angles `latitude`, `longitude`, inclination, and declination are in degrees
  at function boundaries; trigonometric calls convert with `pi/180` in
  `src/MaGICS/geometry.c:50`, `src/MaGICS/geometry.c:73-76`, and
  `src/MaGICS/init.c:136-141`.
- `chi`, `v`, `eta`, `P`, and local probability `proba` are dimensionless.

## 5. Physical constants

Declared constants are in `src/MaGICS/constants.h:1-10`:

| Constant | Value in code | Unit/comment |
| --- | --- | --- |
| `pi` | `3.141592654` | dimensionless |
| `c` | `299792458` | m/s |
| `hbar` | `1.05457E-34` | code comment says "joules per second"; dimensionally this is J s |
| `e` | `1.6021892E-19` | coulomb; also used as J/eV conversion factor |
| `compton_length` | `3.861592642E-13` | meters; not used in the conversion path found here |
| `electron_mass` | `9.109534E-31` | kg |
| `alpha` | `1/137.0` | fine-structure constant |
| `earth_radius` | `6.378140E6` | meters |
| `atmospheric_thickness` | `1E5` | meters |

Derived constant:

- `bcrit = sqr(electron_mass * c) / (e * hbar)` at
  `src/MaGICS/init.c:169-170`. This is the Schwinger critical magnetic field in
  tesla in the code's SI convention.

Other numeric constants affecting conversion:

- Pair conversion step sizes: `dt = 1E-5 s` above `0.1 * earth_radius`, and
  `dt = 5E-7 s` below that altitude in `src/MaGICS/functions.c:225-228` and
  the propagators.
- Low-energy gamma cutoff for immediate atmospheric injection:
  `energy < 5E8` GeV in `src/MaGICS/1D_klepikov_v2.c:14-17` and analogues.
- Atmospheric boundary in propagators often uses `atmospheric_thickness`
  (`1E5 m`), while `conversion_probability()` uses AIRES
  `atm_injection_altitude`.

## 6. Analytical approximations and numerical tables

Local gamma conversion:

- `pair_production_probability()` is explicitly labeled "Erber's analytical
  aproximation" in `src/MaGICS/functions.c:127`.
- Formula in code:
  `1.234E18 * deltatime * (m_e c^2 / (Energy * 1E9 * e)) * K_{1/3}(2/(3 chi))^2`
  at `src/MaGICS/functions.c:131-136`.
- `dbskr3_()` supplies Bessel `K_{nu/3}` values; it comes from the Numerical
  Recipes replacement in `src/NR/dbskr3_fromNR.f`.

Bessel integral / radiation approximations:

- `int_from_x_to_infty_BesselK5_3()` in `src/MaGICS/functions.c:16-49` uses:
  series expansion for `x < 2`, a cubic polynomial for `2 <= x < 5`, and Bessel
  asymptotic expressions for `x >= 5`.
- `kappa()` in `src/MaGICS/functions.c:52-84` repeats the same three-region
  approximation and returns the un-divided function.
- `bremsstrahlung_landau()` in `src/MaGICS/functions.c:139-145` and
  `bremsstrahlung_erber()` in `src/MaGICS/functions.c:147-152` use those
  functions for photon emission, not directly for the reported conversion
  probability.

Numerical sampling:

- `random_electron_energy_fraction()` creates 1025 bins over `v` and samples a
  cumulative sum by binary search in `src/MaGICS/functions.c:101-125`.
- Electron/positron bremsstrahlung in v2 propagators also builds 1025-bin
  arrays over logarithmic energy slices with `energy_slice = 0.9`, for example
  in `src/MaGICS/1D_klepikov_v2.c:90-124`.

Geomagnetic field:

- `cartesian_magnetic_field()` delegates the actual geomagnetic model to AIRES
  `geomagnetic_()` in `src/MaGICS/geometry.c:67-72`; the model/table is not in
  this repository.
- `frenchi()` declares `geomagneticr_()` at `src/MaGICS/geometry.c:215`, but the
  current implementation does not call it.

Trajectory / altitude:

- Propagation is rectilinear for photons in the probability path.
- Altitude calculation uses `heighfactor()` to avoid loss of precision in
  `sqrt(1+x)-1` for small `x`, with a fifth-order series in
  `src/MaGICS/geometry.c:2-18`.

## 7. Random generator and seeds

The generator used by MaGICS is the AIRES routine `urandom_()`, declared in
`src/MaGICS/prototypes.h:5`. There is no RNG implementation and no explicit
seed handling in this repository.

Random uses in the conversion flow:

- Conversion accept/reject: `urandom_() < proba` in
  `src/MaGICS/1D_klepikov_v2.c:33`, `src/MaGICS/1D_klepikov.c:21`,
  `src/MaGICS/1D_erber_v2.c:34`, `src/MaGICS/1D_erber.c:23`, and
  `src/MaGICS/others.c:58`.
- Conversion time inside a step: `timeprod = mytime + urandom_() * dt` in
  `src/MaGICS/1D_klepikov_v2.c:52`, `src/MaGICS/1D_erber_v2.c:45`, and
  `src/MaGICS/others.c:77`.
- Electron energy fraction sampling: `random_electron_energy_fraction()` uses
  `urandom_()` in `src/MaGICS/functions.c:116` and
  `src/MaGICS/functions.c:124`.

Seeds are passed through AIRES input, not through MaGICS code. The example sets
`RandomSeed 0.2560013287` in `example.inp:45-47`. Because `urandom_()` is an
AIRES symbol linked from `-lAires`, that AIRES seed controls MaGICS random
draws.

## 8. Global variables affecting the calculation

Internal globals from `src/MaGICS/vars.h`:

| Global | Lines | Effect |
| --- | --- | --- |
| `debuglvl` | `vars.h:4`, initialized in `init.c:96` | Controls diagnostic printing; can affect output but not physics. |
| `AllConverted` | `vars.h:5`, set in `init.c:207-208` and `init.c:391-399` | Enables forced repeated trials until conversion. |
| `ThereWasNoConversion` | `vars.h:6`, used in `magics.c:69`, propagators | Controls forced-conversion loop and whether particles are returned to AIRES. |
| `FirstInteractionAlready` | `vars.h:6`, propagators | Ensures conversion altitude is written once to AIRES variable 3. |
| `PrimaryParticle` | `vars.h:7`, set by algorithm suffixes in `init.c` | Selects gamma/electron/positron initial particle ID. |
| `bcrit` | `vars.h:9`, set in `init.c:169-170` | Critical field used by `chi()` and `lastbfield`. |
| `geographical_location`, `sintheta_0`, `costheta_0`, `sinphi_0`, `cosphi_0` | `vars.h:10-14`, set in `init.c:131-141` | Site geometry used by rotations. |
| `site_declination`, `sindec_0`, `cosdec_0` | `vars.h:15`, set in `init.c:135-137` | Defines AIRES/local magnetic north rotation. |
| `x_0` | `vars.h:17`, set in `init.c:143-146` | Geocentric site position. |
| `initial_position` | `vars.h:19`, set in `init.c:196-200` | Initial gamma position used by `add_particle()`. |
| `fdate` | `vars.h:21`, set in `init.c:150` | Date passed into `geomagnetic_()`. |
| `mytime` | `vars.h:22`, used throughout propagators | Simulation clock in seconds. |
| `atm_injection_altitude` | `vars.h:26`, set in `init.c:148` | Stop altitude for probability integration. |
| `altitude_0` | `vars.h:27`, set in `init.c:175` | Initial primary altitude. |
| `total_photon_number`, `hard_photons` | `vars.h:24-25` | Counters only; do not change conversion probability. |

Globals from `src/MaGICS/particle_list.c` and AIRES:

| Global | Lines | Effect |
| --- | --- | --- |
| `first_particle` | `particle_list.c:18` | Head of the active particle list; source for probability clone. |
| `lastbfield` | `particle_list.c:20` | Stores last magnetic field direction and normalized strength; used in 3D electron bending in `others.c:156-168`. |
| `shower_number`, `primary_energy`, `default_injection_position`, `injection_depth`, `ground_altitude`, `ground_depth`, `d_ground_inj`, `shower_axis` | `particle_list.c:22-28`, filled by `speistart_()` in `magics.c:64-67` | Initial event, geometry, energy, and direction. |
| `parstring`, `parstringlen` | `particle_list.c:30-33`, filled by `speigetparsc()` in `init.c:105` | Module arguments select algorithm and forced conversion. |
| `realdata`, `intdata`, `shprimcode`, `paquesirve` | `particle_list.c:31-32`, filled by `croinputdata0_()` in `init.c:128` | Site location, field metadata, injection altitude, and date. |
| `irc` | `particle_list.c:43` and `init.c:12` | Return/status code used in AIRES calls. |

AIRES/global input parameters read by name:

- `PEMF_Deblevel`, read with default `0` in `src/MaGICS/init.c:96`.
- `PEMF_init_altit`, read with default `5 * earth_radius` in
  `src/MaGICS/init.c:174-175`.
- Module arguments from `AddSpecialParticle`, parsed in
  `src/MaGICS/init.c:105-122`, select algorithm and `ForceConversion` /
  `NoForceConversion`.

Important AIRES `realdata` indices:

- `realdata[9]`: atmospheric injection altitude, assigned to
  `atm_injection_altitude` in `src/MaGICS/init.c:148`.
- `realdata[20]`, `realdata[21]`: latitude and longitude, assigned in
  `src/MaGICS/init.c:131-134`.
- `realdata[24]`: site declination, assigned in `src/MaGICS/init.c:135`.
- `realdata[26]`: date for geomagnetic model, assigned to `fdate` in
  `src/MaGICS/init.c:150`.
- `realdata[22]`, `realdata[23]`, `realdata[24]` are printed as field,
  inclination, and declination for debug in `src/MaGICS/init.c:155-156`, but
  the actual field in `chi()` is obtained via `geomagnetic_()`.

## 9. Algorithm selection and observed dispatch

`inic()` selects the function pointer `propagate_particle` from the first
module argument:

- No argument defaults to `propagate_1D_klepikov_v2` at
  `src/MaGICS/init.c:221-226`.
- `1Derberaprox` maps to `propagate_1D_erberaprox` at
  `src/MaGICS/init.c:230-252`.
- `1Dklepikov` maps to `propagate_1D_klepikov` at
  `src/MaGICS/init.c:276-298`.
- `3Dklepikov_v2` maps to `propagate_3D_klepikov_v2` at
  `src/MaGICS/init.c:345-367`.
- `UnconvertedPhoton` maps to `propagate_test` at
  `src/MaGICS/init.c:368-374`.

Observed dispatch details that affect runtime behavior:

- `1Derberaprox_v2` currently assigns `propagate_1D_erberaprox`, not
  `propagate_1D_erberaprox_v2`, at `src/MaGICS/init.c:253-258`.
- `1Dklepikov_v2` currently assigns `propagate_1D_erberaprox`, not
  `propagate_1D_klepikov_v2`, at `src/MaGICS/init.c:299-304`.
- The default no-argument path does use `propagate_1D_klepikov_v2`.

## 10. Returned AIRES variables

The README documents return slots in `read.me:97-102`:

- Variable 1: conversion probability, written in `src/MaGICS/magics.c:105-112`.
- Variable 2: number of trials until conversion, written in
  `src/MaGICS/magics.c:168-174`.
- Variable 3: conversion altitude in meters, written when the first conversion
  is sampled, for example in `src/MaGICS/1D_klepikov_v2.c:36-43`.

