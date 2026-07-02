# Pair Prefactor Provenance

This note audits the origin of the MaGICS literal:

```c
1.234E18
```

used in `pair_production_probability()`. No production code is changed.

## A. Evidencia Directa

The literal appears in the active MaGICS implementation at `src/MaGICS/functions.c:134-143`:

```c
result=1.234E18*deltatime*(electron_mass*c*c/(Energy*1E9*e))*
             sqr(dbskr3_(&besselarg,&mu));
```

The only local comment attached to this formula is:

```c
/*Erber's analytical aproximation*/
```

The same literal appears in sibling source trees:

| Path | Evidence |
| --- | --- |
| `MaGICS.1-4-2/src/MaGICS/functions.c:127-136` | same `pair_production_probability()` body and `1.234E18` |
| `MaGICStest/src/MaGICS/functions.c:127-136` | same `pair_production_probability()` body and `1.234E18` |

No local source comment gives a derivation for the number.

Using MaGICS constants:

```text
alpha * m_e c^2 / hbar = 5.666848305750188e18 s^-1
1.234e18 s^-1 = 0.21775772588580714 * alpha * m_e c^2 / hbar
```

Erber 1966 Eq. (3.4), converted to a time rate, gives:

```text
0.16 * alpha * m_e c^2 / hbar
```

Therefore:

```text
C_M / C_Erber = 0.21775772588580714 / 0.16
              = 1.3609857867862947
```

## B. Búsqueda Histórica

Git history was checked with:

```bash
git log --all -S'1.234E18' --oneline --decorate -- src/MaGICS/functions.c
git log --all -S'1.234' --oneline --decorate -- .
git blame -L 134,143 -- src/MaGICS/functions.c
```

Findings:

- `1.234E18` is present in the initial imported file commit `51c85a2` (`version 1.5.0 enviada a MT (files added)`, 2022-09-14).
- `git blame` assigns the literal line to that initial import, so git does not expose the older origin.
- Later commits only add diagnostic instrumentation and do not explain or alter the prefactor.

Sibling directory search under `/home/daniel/2_AREAS/DESARROLLO/MaGICS_TdF` found the same literal in `MaGICS.1-4-2` and `MaGICStest`, but no older explanatory comments.

AIRES was searched in read-only mode. Matches for `0.16`, `0.218`, and `1.234` in AIRES are unrelated shower-model constants, event data, IGRF coefficients, or large tabulated model data. No AIRES formula matching MaGICS magnetic pair production was found.

## C. Candidatos Algebraicos

`tools/analyze_pair_prefactor.py` computes:

```text
C_M = 1.234e18 / (alpha m_e c^2 / hbar)
```

and compares it to simple candidates containing `sqrt(3)`, `pi`, rational factors, and `0.16` scalings. The generated table is `diagnostics/pair_prefactor_candidates.csv`.

Closest simple candidates from the generated search are numerical coincidences only; none has supporting formula evidence with the same `chi`, Bessel argument, units, and polarization convention.

Important rejected/insufficient examples:

| Candidate | Value | Relative difference vs `C_M` | Status |
| --- | ---: | ---: | --- |
| `sqrt(3)/8` | `0.21650635094610965` | about `5.75e-3` | close but not exact; no source formula found |
| `0.16 * (C_M/0.16)` | `0.21775772588580714` | exact by construction | tautology, not provenance |
| `0.16 * sqrt(3)` | `0.27712812921102037` | large mismatch | rejected |
| `3*sqrt(3)/(8*pi)` | `0.20674833578317203` | about `5.05e-2` | not exact; no source formula found |
| `sqrt(3)/(2*pi)` | `0.27566444771089604` | large mismatch | rejected |

No tested simple combination of `sqrt(3)`, `pi`, factors of 2, and low-order rationals reproduced `C_M` exactly with independent formula evidence.

## D. Polarización

The local primary source `erber1966.pdf` presents the approximate total attenuation coefficient in Sec. 3B:

```text
x = (1/2) * (h nu / m c^2) * (H / Hcr)
alpha(x) = (1/2) * (alpha / lambda_c) * (H / Hcr) * T(x)
T(x) ~= 0.16 * x^-1 * K_{1/3}(2/(3x))^2
```

After converting to a time rate:

```text
Gamma_E = 0.16 * alpha * m_e c^2 / hbar
          * (m_e c^2 / E)
          * K_{1/3}(2/(3x))^2
```

The page does not provide a polarization-resolved split for this approximate expression. No local source was found that changes the coefficient from `0.16` to `0.21775772588580714` through a documented polarization average or polarization-specific rate while preserving:

- `chi = x`;
- Bessel argument `2/(3 chi)`;
- a time-rate convention;
- the same `K_{1/3}^2` approximation.

Therefore a polarization origin for `1.234E18` is not demonstrated.

## E. Fuentes Bibliográficas

Local sources checked:

- `erber1966.pdf`: primary source; supports `0.16`, not `0.21775772588580714`.
- `2005ICRC....9....1B.pdf`: local PDF exists, but `pdftotext` extracts only the ADS identifier (`2005ICRC....9....1B`) and no readable formula text.
- `docs/pair_conversion_physics_audit.md` and `docs/erber_formula_comparison.md`: local audit documents; no independent provenance for `1.234E18`.
- AIRES manuals and source in `/home/daniel/2_AREAS/DESARROLLO/MaGICS_TdF/aires` and `/home/daniel/2_AREAS/DESARROLLO/MaGICS_TdF/19-04-08`: no matching magnetic-pair prefactor formula found.
- MaGICS sibling trees `MaGICS.1-4-2` and `MaGICStest`: same literal, no derivation.

No local Klepikov or PRESHOWER source formula was found that provides the `1.234E18` coefficient with the MaGICS `chi` and Bessel argument.

## F. Hipótesis Descartadas

### Direct Erber Eq. (3.4)

Rejected as exact origin of `1.234E18`.

Erber gives coefficient `0.16`; MaGICS implies `0.21775772588580714`. Since `chi` and the Bessel argument match, this is a pure prefactor discrepancy:

```text
1.234e18 / [0.16 * alpha * m_e c^2 / hbar] = 1.3609857867862947
```

### Simple `sqrt(3)` / `pi` factor

No exact simple factor was found. The nearest obvious candidate, `sqrt(3)/8`, is close but differs by about `0.575%` and lacks documentary support.

### Unit conversion artifact

Rejected. Unit conversion issues in `chi()` cancel GeV/nT factors and do not affect the standalone rate prefactor once written as `C_M * alpha m_e c^2/hbar`.

### Temporal vs length-rate conversion

Rejected as an explanation. Converting Erber's length attenuation to a time rate multiplies by `c` and produces the scale `alpha m_e c^2/hbar`; the remaining coefficient is still `0.16`.

### AIRES source

Rejected. AIRES contains unrelated `0.16` and `1.234...` numerical data but no local magnetic-pair production prefactor formula matching MaGICS.

## G. Conclusión

The origin of `1.234E18` is not resolved from local evidence.

What is demonstrated:

- The literal is direct code data in MaGICS and predates the current git history as an imported line.
- It is also present in sibling MaGICS trees.
- It does not equal Erber 1966 Eq. (3.4)'s coefficient under the already demonstrated shared `chi` and Bessel-argument convention.
- It corresponds to `C_M = 0.21775772588580714` in units of `alpha m_e c^2/hbar`.

What is not demonstrated:

- A published formula with coefficient `0.21775772588580714`.
- A polarization-resolved normalization yielding this coefficient.
- A simple exact algebraic combination of `sqrt(3)`, `pi`, powers of two, and Erber's `0.16`.
- A historical commit or note explaining why this value was chosen.

The most defensible classification is:

```text
origin: no resuelto
status vs Erber 1966 Eq. (3.4): demonstrated prefactor mismatch
bug status: not proven
```

It could be a transcription error, an undocumented alternate normalization, or a coefficient copied from a source not present locally. The available evidence does not distinguish these.

## H. Nivel de Confianza

| Claim | Confidence | Basis |
| --- | --- | --- |
| `1.234E18` is hard-coded in MaGICS | High | direct code and sibling trees |
| Git history does not expose older derivation | High | `git log -S`, `git blame` |
| Erber 1966 Eq. (3.4) gives `0.16`, not `0.2177577` | High | local primary PDF |
| Unit conversion does not explain the prefactor | High | algebraic separation of `chi` and rate scale |
| Polarization explains the coefficient | Low | no local supporting formula |
| Transcription error | Plausible but unproven | no direct evidence |
| Undocumented alternate source | Plausible but unproven | no local source found |
