# Legacy MaGICS Inventory

This inventory was prepared before moving any legacy source. The MaGICS code
remains in `src/MaGICS` and `src/NR`; `legacy/` is reserved for future staged
migration but is intentionally empty in this bootstrap.

| Component | Legacy evidence | Classification | Rationale |
| --- | --- | --- | --- |
| Geometry | `src/MaGICS/geometry.c`, `docs/geomagnetic_interface.md` | Reimplement with the same formula | The coordinate transforms, altitude calculation and `B_perp` projection are conceptually useful, but the code is tied to AIRES globals, hidden units and mutable `lastbfield`. |
| AIRES interface | `src/MaGICS/magics.c`, `src/MaGICS/init.c`, `src/MaGICS/prototypes.h` | Reuse only as reference | The call sequence, returned variables and special-primary lifecycle are valuable documentation. The new core must not include AIRES headers or symbols. |
| Pair production | `src/MaGICS/functions.c`, `docs/erber_formula_comparison.md`, `docs/pair_prefactor_provenance.md` | Reimplement with the same source formula | The Bessel argument and Erber structure are useful, but the literal `1.234E18` is undocumented. The new bootstrap uses Erber 1966 coefficient `0.16` built from constants. |
| Pair energy distribution | `pair_production_spectrum()` and `random_electron_energy_fraction()` in `src/MaGICS/functions.c` | Reuse only as reference | The Klepikov spectrum and sampling strategy are important, but the current implementation uses globals, AIRES RNG and ad hoc bins. |
| Magnetic emission | `bremsstrahlung_erber()`, `bremsstrahlung_landau()`, propagators in `src/MaGICS/1D_*` | Reuse only as reference | Formula choices and historical behavior should be preserved for comparison, but transport and emission need a clearer optical-depth design. |
| Particle management | `src/MaGICS/particle_list.c` | Discard | Intrusive linked-list state and raw allocation are unsuitable for a reusable core. Keep only the concept of explicit particle records. |
| IGRF / geomagnetic model | AIRES `geomagnetic_()` usage documented in `docs/geomagnetic_interface.md` | Reuse only as reference | The implementation is outside this repository. Future standalone mode needs its own field-provider interface and eventually an explicit IGRF source or dependency. |
| Tests and historical cases | `tests/`, `diagnostics/`, `docs/*audit*.md`, `erber1966.pdf` | Reuse almost directly | These are already independent evidence and should remain regression/reference material. Some Python tests depend on external diagnostic CSV paths and need later cleanup. |

## Notable Legacy Risks

- AIRES symbols appear throughout the legacy executable path.
- Physical units are implicit; `chi()` relies on a GeV/nT numerical
  cancellation.
- The RNG is AIRES-provided through `urandom_()`.
- There is mutable global state in `vars.h` and `particle_list.c`.
- Some algorithm dispatch paths in `init.c` appear inconsistent with their
  printed names.
- The Erber prefactor literal `1.234E18` is not derived locally.
