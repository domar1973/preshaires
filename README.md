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
- a minimal CLI command;
- unit tests for the bootstrap formula and unit conversions;
- an empty, compilable AIRES adapter placeholder.

It does not implement transport, cascades, magnetic emission, pair energy
sampling, Python bindings or real AIRES coupling.

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

## Historical MaGICS

The legacy MaGICS_TdF code remains under `src/MaGICS` and `src/NR`, with the
original `read.me`, `MakeMaGICS`, diagnostics and audit notes preserved. New
code should not include AIRES headers or symbols outside `adapters/aires`.

## Licensing

No explicit license file was found during bootstrap. See
`docs/architecture/licensing.md` before redistributing or relicensing any part
of this repository.
