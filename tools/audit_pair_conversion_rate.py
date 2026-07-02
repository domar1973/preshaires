#!/usr/bin/env python3
"""Independent audit of MaGICS magnetic pair-conversion step rates."""

import argparse
import math
import os
from pathlib import Path

os.makedirs("/tmp/matplotlib", exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import kv


ELEMENTARY_CHARGE = 1.6021892e-19
ELECTRON_MASS = 9.109534e-31
C_LIGHT = 299792458.0
HBAR = 1.05457e-34
ALPHA = 1.0 / 137.0
ERBER_PAIR_COEFFICIENT = 0.16
RATE_PREFACTOR = ERBER_PAIR_COEFFICIENT * ALPHA * ELECTRON_MASS * C_LIGHT**2 / HBAR


def bessel_argument(chi):
    chi = np.asarray(chi, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 / (3.0 * chi)


def gamma_rate(energy_gev, chi):
    energy_gev = np.asarray(energy_gev, dtype=float)
    chi = np.asarray(chi, dtype=float)
    e_j = energy_gev * 1.0e9 * ELEMENTARY_CHARGE
    arg = bessel_argument(chi)
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        kval = kv(1.0 / 3.0, arg)
        return RATE_PREFACTOR * (ELECTRON_MASS * C_LIGHT**2 / e_j) * kval**2


def relative_error(reference, observed):
    reference = np.asarray(reference, dtype=float)
    observed = np.asarray(observed, dtype=float)
    denominator = np.maximum(np.abs(reference), np.finfo(float).tiny)
    return np.abs(observed - reference) / denominator


def finite_stats(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {
            "max": math.nan,
            "mean": math.nan,
            "median": math.nan,
            "rms": math.nan,
            "p50": math.nan,
            "p90": math.nan,
            "p99": math.nan,
            "p100": math.nan,
        }
    return {
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "rms": float(np.sqrt(np.mean(finite**2))),
        "p50": float(np.percentile(finite, 50)),
        "p90": float(np.percentile(finite, 90)),
        "p99": float(np.percentile(finite, 99)),
        "p100": float(np.percentile(finite, 100)),
    }


def analyze_pair_rate(csv_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_path)

    required = [
        "energy_GeV",
        "chi",
        "dt_s",
        "step_probability",
        "accumulated_survival",
        "accumulated_probability",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required CSV columns: {missing}")

    rate = gamma_rate(df["energy_GeV"].to_numpy(), df["chi"].to_numpy())
    reference_step = rate * df["dt_s"].to_numpy()
    abs_diff = np.abs(reference_step - df["step_probability"].to_numpy())
    rel_diff = relative_error(reference_step, df["step_probability"].to_numpy())
    arg = bessel_argument(df["chi"].to_numpy())

    comparison = df.copy()
    comparison["bessel_argument"] = arg
    comparison["python_gamma_rate_s^-1"] = rate
    comparison["reference_step_probability"] = reference_step
    comparison["absolute_difference"] = abs_diff
    comparison["relative_difference"] = rel_diff
    comparison.to_csv(output_dir / "pair_rate_row_comparison.csv", index=False)

    finite_abs = abs_diff[np.isfinite(abs_diff)]
    finite_rel = rel_diff[np.isfinite(rel_diff)]
    max_abs_idx = int(np.nanargmax(np.where(np.isfinite(abs_diff), abs_diff, -np.inf)))
    max_rel_idx = int(np.nanargmax(np.where(np.isfinite(rel_diff), rel_diff, -np.inf)))

    step = df["step_probability"].to_numpy()
    nonphysical_mask = (
        (df["chi"].to_numpy() <= 0)
        | (df["dt_s"].to_numpy() <= 0)
        | (df["energy_GeV"].to_numpy() <= 0)
        | (step < 0)
    )
    nan_count = int(
        df[required].isna().sum().sum()
        + np.isnan(rate).sum()
        + np.isnan(reference_step).sum()
        + np.isnan(rel_diff).sum()
    )
    inf_count = int(
        np.isinf(df[required].to_numpy(dtype=float)).sum()
        + np.isinf(rate).sum()
        + np.isinf(reference_step).sum()
        + np.isinf(rel_diff).sum()
    )

    p_product = float(1.0 - np.prod(1.0 - step))
    p_exp_linear = float(1.0 - np.exp(-np.sum(step)))
    tau_reference = float(np.sum(reference_step))
    p_reference = float(1.0 - np.exp(-tau_reference))
    p_csv_final = float(df["accumulated_probability"].iloc[-1])

    abs_stats = finite_stats(abs_diff)
    rel_stats = finite_stats(rel_diff)

    def rel_to_csv(value):
        denom = max(abs(p_csv_final), np.finfo(float).tiny)
        return abs(value - p_csv_final) / denom

    summary_path = output_dir / "pair_rate_summary.md"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Pair Conversion Rate Audit\n\n")
        handle.write(f"Input CSV: `{csv_path}`\n\n")
        handle.write(f"Rows: {len(df)}\n\n")
        handle.write("## Row-wise Step Probability Comparison\n\n")
        handle.write(f"- Maximum absolute error: `{abs_stats['max']:.17e}`\n")
        handle.write(f"- Maximum relative error: `{rel_stats['max']:.17e}`\n")
        handle.write(f"- Mean relative error: `{rel_stats['mean']:.17e}`\n")
        handle.write(f"- Median relative error: `{rel_stats['median']:.17e}`\n")
        handle.write(f"- RMS relative error: `{rel_stats['rms']:.17e}`\n")
        handle.write(f"- Relative-error p50: `{rel_stats['p50']:.17e}`\n")
        handle.write(f"- Relative-error p90: `{rel_stats['p90']:.17e}`\n")
        handle.write(f"- Relative-error p99: `{rel_stats['p99']:.17e}`\n")
        handle.write(f"- Relative-error p100: `{rel_stats['p100']:.17e}`\n")
        handle.write(f"- NaN count across audited values: `{nan_count}`\n")
        handle.write(f"- Inf count across audited values: `{inf_count}`\n")
        handle.write(f"- Non-physical row count: `{int(nonphysical_mask.sum())}`\n\n")
        handle.write("## Maximum Error Rows\n\n")
        for label, idx in [("absolute", max_abs_idx), ("relative", max_rel_idx)]:
            row = comparison.iloc[idx]
            handle.write(f"### Maximum {label} error\n\n")
            handle.write(f"- Row index: `{idx}`\n")
            handle.write(f"- `energy_GeV`: `{row['energy_GeV']:.17e}`\n")
            handle.write(f"- `chi`: `{row['chi']:.17e}`\n")
            handle.write(f"- `dt_s`: `{row['dt_s']:.17e}`\n")
            handle.write(f"- Bessel argument `2/(3 chi)`: `{row['bessel_argument']:.17e}`\n")
            handle.write(f"- MaGICS `step_probability`: `{row['step_probability']:.17e}`\n")
            handle.write(
                f"- Python `reference_step_probability`: "
                f"`{row['reference_step_probability']:.17e}`\n"
            )
            handle.write(f"- Absolute difference: `{row['absolute_difference']:.17e}`\n")
            handle.write(f"- Relative difference: `{row['relative_difference']:.17e}`\n\n")
        handle.write("## Integrated Probability Reconstruction\n\n")
        handle.write(f"- CSV final `accumulated_probability`: `{p_csv_final:.17e}`\n")
        handle.write(f"- `P_product = 1 - product(1 - p_i)`: `{p_product:.17e}`\n")
        handle.write(f"- `P_exp_linear = 1 - exp(-sum(p_i))`: `{p_exp_linear:.17e}`\n")
        handle.write(f"- `tau_reference = sum(gamma_rate_i * dt_i)`: `{tau_reference:.17e}`\n")
        handle.write(f"- `P_reference = 1 - exp(-tau_reference)`: `{p_reference:.17e}`\n\n")
        handle.write("| Quantity | Absolute difference vs CSV final | Relative difference vs CSV final |\n")
        handle.write("| --- | ---: | ---: |\n")
        for label, value in [
            ("P_product", p_product),
            ("P_exp_linear", p_exp_linear),
            ("P_reference", p_reference),
        ]:
            handle.write(
                f"| `{label}` | `{abs(value - p_csv_final):.17e}` | "
                f"`{rel_to_csv(value):.17e}` |\n"
            )
        handle.write("\n")
        handle.write("## Interpretation\n\n")
        handle.write(
            "`P_product` should reproduce the CSV final accumulated probability "
            "up to floating-point roundoff if the CSV contains every step used by "
            "MaGICS. Differences between `P_product` and `P_exp_linear` measure "
            "the effect of replacing the discrete product of linear probabilities "
            "with an exponential survival model. Differences involving "
            "`P_reference` also include any row-wise disagreement between the "
            "independent SciPy Bessel evaluation and the MaGICS step probability.\n\n"
        )
        product_exp_gap = abs(p_product - p_exp_linear) / max(abs(p_product), np.finfo(float).tiny)
        handle.write(
            f"The relative gap between `P_product` and `P_exp_linear` is "
            f"`{product_exp_gap:.17e}`. A discrepancy of order 10-20 percent "
            "would require this gap, or the independent rate gap, to be of that "
            "same order.\n"
        )

    plt.figure()
    finite_plot = np.isfinite(rel_diff)
    plt.semilogy(np.arange(len(rel_diff))[finite_plot], rel_diff[finite_plot], ".")
    plt.xlabel("CSV row")
    plt.ylabel("relative error")
    plt.tight_layout()
    plt.savefig(output_dir / "pair_rate_relative_error.png", dpi=160)
    plt.close()

    plt.figure()
    chi = df["chi"].to_numpy()
    magics_rate = step / df["dt_s"].to_numpy()
    mask = (chi > 0) & np.isfinite(chi) & np.isfinite(rate) & np.isfinite(magics_rate)
    positive_rates = (rate > 0) & (magics_rate > 0)
    if np.any(mask & positive_rates):
        plt.loglog(chi[mask & positive_rates], magics_rate[mask & positive_rates], ".", label="MaGICS")
        plt.loglog(chi[mask & positive_rates], rate[mask & positive_rates], ".", label="Python")
    else:
        plt.plot(chi[mask], magics_rate[mask], ".", label="MaGICS")
        plt.plot(chi[mask], rate[mask], ".", label="Python")
    plt.xlabel("chi")
    plt.ylabel("rate s^-1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "pair_rate_vs_chi.png", dpi=160)
    plt.close()

    return {
        "rows": len(df),
        "max_abs_error": abs_stats["max"],
        "max_rel_error": rel_stats["max"],
        "mean_rel_error": rel_stats["mean"],
        "median_rel_error": rel_stats["median"],
        "rms_rel_error": rel_stats["rms"],
        "p_csv_final": p_csv_final,
        "p_product": p_product,
        "p_exp_linear": p_exp_linear,
        "p_reference": p_reference,
        "tau_reference": tau_reference,
        "summary_path": str(summary_path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="MaGICS diagnostic CSV")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("diagnostics"),
        help="Directory for generated audit artifacts",
    )
    args = parser.parse_args()
    result = analyze_pair_rate(args.csv, args.output_dir)
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
