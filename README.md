# preshaires

preshaires is the architectural bootstrap for a geomagnetic preshower simulator
derived from the historical MaGICS_TdF repository.

The goal is a shared physical core that can run as a standalone library/CLI and
later through an AIRES adapter. The current repository still contains the
MaGICS legacy code in place; it has not been moved or deleted.

## Status

Experimental. This stage implements only:

- explicit energy, magnetic-field, length, time and rate wrappers;
- photon `chi` for a transverse magnetic field;
- the local Erber 1966 approximate pair-production rate;
- deterministic conversion probability over a prescribed straight trajectory;
- minimal CLI commands for local rate and trajectory probability;
- unit tests for the bootstrap formula and unit conversions;
- an empty, compilable AIRES adapter placeholder.

It does not implement Monte Carlo interaction sampling, cascades, magnetic
emission, pair energy sampling, Earth geometry, IGRF, Python bindings or real
AIRES coupling.

## Build

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## CLI

```sh
./build/preshaires rate --energy-eV 7e19 --Bperp-T 2.115138e-5
```

The command prints the model name, photon `chi` and local rate in `s^-1`.

For a prescribed straight path, preshaires integrates optical depth:

```text
tau = integral_0^L [Gamma(E, B_perp(s)) / c] ds
P = 1 - exp(-tau)
```

`Gamma` is the local Erber rate per unit time, so division by `c` converts it
to attenuation per unit length before integrating over path length.

Uniform field example:

```sh
./build/preshaires probability \
  --energy-eV 7e19 \
  --length-m 1e7 \
  --direction 0,0,-1 \
  --uniform-field-T 2.115138e-5,0,0
```

Tabulated field example:

```sh
./build/preshaires probability \
  --energy-eV 7e19 \
  --length-m 1e7 \
  --direction 0,0,-1 \
  --field-table tests/reference/uniform_field.csv
```

The table format is:

```text
s_m,Bx_T,By_T,Bz_T
0,2.115138e-5,0,0
10000000,2.115138e-5,0,0
```

The probability command is deterministic. It does not sample a conversion
point and does not create an electron-positron pair.

## Historical MaGICS

The legacy MaGICS_TdF code remains under `src/MaGICS` and `src/NR`, with the
original `read.me`, `MakeMaGICS`, diagnostics and audit notes preserved. New
code should not include AIRES headers or symbols outside `adapters/aires`.

## Licensing

No explicit license file was found during bootstrap. See
`docs/architecture/licensing.md` before redistributing or relicensing any part
of this repository.
