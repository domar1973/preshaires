# Erber Convention Numerical Comparison

Reference evaluated: Erber 1966 Eq. (3.4) using the analytic approximation Eq. (3.3d), with the photon parameter `x = chi_M`.

- MaGICS coefficient `C_M`: `2.17757725885807163e-01`
- Erber coefficient `C_E`: `1.60000000000000003e-01`
- Prefactor ratio `C_M / C_E`: `1.36098578678629467e+00`
- Finite numerical ratio range: [`1.36098578678629445e+00`, `1.36098578678629512e+00`]

## Diagnostic Case

- `E_eV`: `7.00000000000000000e+19`
- `Bperp_T`: `2.11513800000000008e-05`
- `chi_M`: `3.28199918892518316e-01`
- `chi_source`: `3.28199918892518316e-01`
- `rate_MaGICS`: `1.13251095047984478e+02`
- `rate_source`: `8.32125479542332727e+01`
- `ratio`: `1.36098578678629512e+00`
- `relative_difference`: `3.60985786786295060e-01`

The ratio is constant wherever both rates are finite because the Bessel argument and energy dependence are algebraically identical; only the dimensionless prefactor differs.
