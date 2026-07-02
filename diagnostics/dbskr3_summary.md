# `dbskr3_(x, 1)` Audit

The local Fortran implementation sets `nu = abs(nu3) / 3.d0`, calls `bessik(x, nu, ri, rk, rip, rkp)`, and returns `rk`. For `nu3 = 1`, this is `K_{1/3}(x)` as implemented by the Numerical Recipes `bessik` routine.

Rows: `600`

- Maximum relative error: `6.15114933774050510e+00`
- Median relative error: `9.97783246270804004e-08`
- 99th percentile relative error: `4.02799480150904925e-03`
- Maximum absolute error: `5.41665214031955555e-04`
- Non-finite Fortran values: `0`
- Non-finite SciPy values: `0`

## Maximum Relative Error Row

- Row index: `117`
- `chi`: `9.47605805704949825e-04`
- `x`: `7.03527418946863804e+02`
- `dbskr3_x_1`: `1.36867615911802821e-307`
- `scipy_kv_1_3`: `0.00000000000000000e+00`
- Relative difference: `6.15114933774050510e+00`

Relative error above `1e-8` appears for `x` in [`2.00157182161985325e+00`, `7.17180171099230506e+02`].
Fortran underflow to exact zero appears for `x >= 7.45285658699168380e+02`.

## Diagnostic CSV Range

- CSV `chi` range: [`2.29583822493718032e-05`, `3.28704149406820967e-01`]
- CSV `x = 2/(3 chi)` range: [`2.02816626400893441e+00`, `2.90380506529334561e+04`]
- Probe rows overlapping CSV `x` range: `422` of `600`
- Maximum relative error on overlapping probe rows: `6.15114933774050510e+00`
- Maximum absolute error on overlapping probe rows: `5.01593388740020907e-04`
The CSV reaches below the requested probe grid lower bound (`1.00000000000000005e-04`), where both SciPy and Fortran rates are already effectively zero for the conversion-rate audit.
