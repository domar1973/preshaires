# Pair Conversion Rate Audit

Input CSV: `/home/daniel/2_AREAS/DESARROLLO/MaGICS_TdF/TEST/docs/conversion_probability_diagnostics.csv`

Rows: 40532

## Row-wise Step Probability Comparison

- Maximum absolute error: `2.64897956966162937e-04`
- Maximum relative error: `3.73440983685451866e-01`
- Mean relative error: `1.45355441837933169e-01`
- Median relative error: `0.00000000000000000e+00`
- RMS relative error: `2.30896670427605893e-01`
- Relative-error p50: `0.00000000000000000e+00`
- Relative-error p90: `3.73255685372911816e-01`
- Relative-error p99: `3.73440056443386070e-01`
- Relative-error p100: `3.73440983685451866e-01`
- NaN count across audited values: `0`
- Inf count across audited values: `0`
- Non-physical row count: `0`

## Maximum Error Rows

### Maximum absolute error

- Row index: `34047`
- `energy_GeV`: `7.00000000000000000e+10`
- `chi`: `3.18144714679375662e-01`
- `dt_s`: `1.00000000000000008e-05`
- Bessel argument `2/(3 chi)`: `2.09548245155833968e+00`
- MaGICS `step_probability`: `9.75382093344335848e-04`
- Python `reference_step_probability`: `7.10484136378172911e-04`
- Absolute difference: `2.64897956966162937e-04`
- Relative difference: `3.72841480059681984e-01`

### Maximum relative error

- Row index: `39409`
- `energy_GeV`: `7.00000000000000000e+10`
- `chi`: `3.28704149406820967e-01`
- `dt_s`: `5.00000000000000083e-07`
- Bessel argument `2/(3 chi)`: `2.02816626400893441e+00`
- MaGICS `step_probability`: `5.75850761958482540e-05`
- Python `reference_step_probability`: `4.19275941812411357e-05`
- Absolute difference: `1.56574820146071183e-05`
- Relative difference: `3.73440983685451866e-01`

## Integrated Probability Reconstruction

- CSV final `accumulated_probability`: `5.25407086018462310e-01`
- `P_product = 1 - product(1 - p_i)`: `5.25407086018462088e-01`
- `P_exp_linear = 1 - exp(-sum(p_i))`: `5.25347205266517792e-01`
- `tau_reference = sum(gamma_rate_i * dt_i)`: `5.43086753163127822e-01`
- `P_reference = 1 - exp(-tau_reference)`: `4.19047774265825512e-01`

| Quantity | Absolute difference vs CSV final | Relative difference vs CSV final |
| --- | ---: | ---: |
| `P_product` | `2.22044604925031308e-16` | `4.22614408586847416e-16` |
| `P_exp_linear` | `5.98807519445188063e-05` | `1.13970202416369115e-04` |
| `P_reference` | `1.06359311752636798e-01` | `2.02432198923368584e-01` |

## Interpretation

`P_product` should reproduce the CSV final accumulated probability up to floating-point roundoff if the CSV contains every step used by MaGICS. Differences between `P_product` and `P_exp_linear` measure the effect of replacing the discrete product of linear probabilities with an exponential survival model. Differences involving `P_reference` also include any row-wise disagreement between the independent SciPy Bessel evaluation and the MaGICS step probability.

The relative gap between `P_product` and `P_exp_linear` is `1.13970202415946548e-04`. A discrepancy of order 10-20 percent would require this gap, or the independent rate gap, to be of that same order.
