# Known Issues

## Resolved: Erber pair-production prefactor

Resolved on branch `fix/erber-pair-prefactor`. The old
`pair_production_probability()` prefactor `1.234E18 s^-1` was a manual
calculation error confirmed by the original author. It has been replaced with
the Erber 1966 Eq. (3.4) time-rate prefactor
`0.16 * alpha * electron_mass * c * c / hbar`. See
`docs/erber_prefactor_bugfix.md`.

## `inic()` dispatch mismatches

These issues are recorded only; they are not corrected in this change.

| Requested algorithm | File and line | Observed dispatch | Expected dispatch |
| --- | --- | --- | --- |
| `1Derberaprox_v2` | `src/MaGICS/init.c:253-255` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Derberaprox_v2+` | `src/MaGICS/init.c:260-263` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Derberaprox_v2-` | `src/MaGICS/init.c:268-271` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Dklepikov_v2` | `src/MaGICS/init.c:299-302` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |
| `1Dklepikov_v2+` | `src/MaGICS/init.c:306-309` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |
| `1Dklepikov_v2-` | `src/MaGICS/init.c:314-317` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |

The default branch for no arguments does select `propagate_1D_klepikov_v2` at `src/MaGICS/init.c:221-225`, so the mismatch affects explicit `*_v2` arguments.

## Implicit `chi()` unit cancellation

`chi()` in `src/MaGICS/geometry.c:196` uses `particle.energy` numerically in GeV and `geomagnetic_()` output numerically in nT. The missing `1e9` energy conversion and the nT-to-T factor cancel, so the current `chi` is numerically correct. This is fragile: changing only one of those unit conventions would change `chi` by `1e9`. The formula is intentionally not changed here.
