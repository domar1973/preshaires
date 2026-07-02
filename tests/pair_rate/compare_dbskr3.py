#!/usr/bin/env python3
"""Compare MaGICS dbskr3(x, 1) probe output against scipy.special.kv."""

import argparse
import os
import subprocess
from pathlib import Path

os.makedirs("/tmp/matplotlib", exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import kv


def relative_error(reference, observed):
    denominator = np.maximum(np.abs(reference), np.finfo(float).tiny)
    return np.abs(observed - reference) / denominator


def run_probe(probe):
    proc = subprocess.run([str(probe)], check=True, text=True, capture_output=True)
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    rows = []
    for line in lines:
        if line.startswith("chi"):
            continue
        parts = [float(part) for part in line.split(",")]
        rows.append(parts)
    return pd.DataFrame(rows, columns=["chi", "x", "dbskr3_x_1"])


def compare(probe, output_dir, csv_path=None, pair_csv=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    if csv_path is None:
        df = run_probe(probe)
    else:
        df = pd.read_csv(csv_path)
    df["scipy_kv_1_3"] = kv(1.0 / 3.0, df["x"].to_numpy())
    df["absolute_difference"] = np.abs(df["dbskr3_x_1"] - df["scipy_kv_1_3"])
    df["relative_difference"] = relative_error(
        df["scipy_kv_1_3"].to_numpy(), df["dbskr3_x_1"].to_numpy()
    )
    out_csv = output_dir / "dbskr3_comparison.csv"
    df.to_csv(out_csv, index=False)

    finite_rel = df["relative_difference"].to_numpy()
    finite_rel = finite_rel[np.isfinite(finite_rel)]
    max_idx = int(
        np.nanargmax(
            np.where(np.isfinite(df["relative_difference"]), df["relative_difference"], -np.inf)
        )
    )
    max_row = df.iloc[max_idx]
    summary = output_dir / "dbskr3_summary.md"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# `dbskr3_(x, 1)` Audit\n\n")
        handle.write("The local Fortran implementation sets `nu = abs(nu3) / 3.d0`, calls `bessik(x, nu, ri, rk, rip, rkp)`, and returns `rk`. For `nu3 = 1`, this is `K_{1/3}(x)` as implemented by the Numerical Recipes `bessik` routine.\n\n")
        handle.write(f"Rows: `{len(df)}`\n\n")
        handle.write(f"- Maximum relative error: `{np.max(finite_rel):.17e}`\n")
        handle.write(f"- Median relative error: `{np.median(finite_rel):.17e}`\n")
        handle.write(f"- 99th percentile relative error: `{np.percentile(finite_rel, 99):.17e}`\n")
        handle.write(f"- Maximum absolute error: `{df['absolute_difference'].max():.17e}`\n")
        handle.write(f"- Non-finite Fortran values: `{int((~np.isfinite(df['dbskr3_x_1'])).sum())}`\n")
        handle.write(f"- Non-finite SciPy values: `{int((~np.isfinite(df['scipy_kv_1_3'])).sum())}`\n\n")
        handle.write("## Maximum Relative Error Row\n\n")
        handle.write(f"- Row index: `{max_idx}`\n")
        handle.write(f"- `chi`: `{max_row['chi']:.17e}`\n")
        handle.write(f"- `x`: `{max_row['x']:.17e}`\n")
        handle.write(f"- `dbskr3_x_1`: `{max_row['dbskr3_x_1']:.17e}`\n")
        handle.write(f"- `scipy_kv_1_3`: `{max_row['scipy_kv_1_3']:.17e}`\n")
        handle.write(f"- Relative difference: `{max_row['relative_difference']:.17e}`\n\n")
        high = df[df["relative_difference"] > 1e-8]
        if high.empty:
            handle.write("No relative-error growth above `1e-8` was observed on the probe grid.\n")
        else:
            handle.write(
                "Relative error above `1e-8` appears for "
                f"`x` in [`{high['x'].min():.17e}`, `{high['x'].max():.17e}`].\n"
            )
        zero_fortran = df[df["dbskr3_x_1"] == 0]
        if zero_fortran.empty:
            handle.write("No Fortran underflow to exact zero was observed on the probe grid.\n")
        else:
            handle.write(
                "Fortran underflow to exact zero appears for "
                f"`x >= {zero_fortran['x'].min():.17e}`.\n"
            )
        if pair_csv is not None:
            pair = pd.read_csv(pair_csv)
            pair_x = 2.0 / (3.0 * pair["chi"].to_numpy())
            handle.write("\n## Diagnostic CSV Range\n\n")
            handle.write(f"- CSV `chi` range: [`{pair['chi'].min():.17e}`, `{pair['chi'].max():.17e}`]\n")
            handle.write(f"- CSV `x = 2/(3 chi)` range: [`{np.nanmin(pair_x):.17e}`, `{np.nanmax(pair_x):.17e}`]\n")
            in_csv_x = df[(df["x"] >= np.nanmin(pair_x)) & (df["x"] <= np.nanmax(pair_x))]
            if in_csv_x.empty:
                handle.write("The probe grid does not overlap the CSV Bessel-argument range.\n")
            else:
                handle.write(
                    "- Probe rows overlapping CSV `x` range: "
                    f"`{len(in_csv_x)}` of `{len(df)}`\n"
                )
                handle.write(
                    "- Maximum relative error on overlapping probe rows: "
                    f"`{in_csv_x['relative_difference'].max():.17e}`\n"
                )
                handle.write(
                    "- Maximum absolute error on overlapping probe rows: "
                    f"`{in_csv_x['absolute_difference'].max():.17e}`\n"
                )
            if pair["chi"].min() < df["chi"].min():
                handle.write(
                    "The CSV reaches below the requested probe grid lower bound "
                    f"(`{df['chi'].min():.17e}`), where both SciPy and Fortran rates "
                    "are already effectively zero for the conversion-rate audit.\n"
                )

    plt.figure()
    mask = np.isfinite(df["relative_difference"])
    plt.loglog(df.loc[mask, "x"], df.loc[mask, "relative_difference"], ".")
    plt.xlabel("x")
    plt.ylabel("relative error")
    plt.tight_layout()
    plt.savefig(output_dir / "dbskr3_relative_error.png", dpi=160)
    plt.close()

    return {
        "rows": len(df),
        "max_relative_error": float(np.max(finite_rel)),
        "summary_path": str(summary),
        "comparison_csv": str(out_csv),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, default=Path("tests/pair_rate/dbskr3_probe"))
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("diagnostics"))
    parser.add_argument("--csv", type=Path, default=None, help="Existing probe CSV to compare")
    parser.add_argument("--pair-csv", type=Path, default=None, help="MaGICS diagnostic CSV")
    args = parser.parse_args()
    result = compare(args.probe, args.output_dir, args.csv, args.pair_csv)
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
