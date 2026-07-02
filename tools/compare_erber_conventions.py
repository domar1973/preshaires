#!/usr/bin/env python3
"""Compare MaGICS magnetic pair-conversion rate with Erber's convention."""

import argparse
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
RATE_PREFACTOR_MAGICS = 1.234e18
COEFFICIENT_ERBER = 0.16
BCRIT_T = (ELECTRON_MASS * C_LIGHT) ** 2 / (ELEMENTARY_CHARGE * HBAR)


def chi_magics_physical(energy_ev, bperp_t):
    energy_j = np.asarray(energy_ev, dtype=float) * ELEMENTARY_CHARGE
    bperp_t = np.asarray(bperp_t, dtype=float)
    return 0.5 * (energy_j / (ELECTRON_MASS * C_LIGHT**2)) * (bperp_t / BCRIT_T)


def magics_rate(energy_ev, bperp_t):
    chi_m = chi_magics_physical(energy_ev, bperp_t)
    energy_j = np.asarray(energy_ev, dtype=float) * ELEMENTARY_CHARGE
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        k = kv(1.0 / 3.0, 2.0 / (3.0 * chi_m))
        return RATE_PREFACTOR_MAGICS * (ELECTRON_MASS * C_LIGHT**2 / energy_j) * k**2


def erber1966_rate(energy_ev, bperp_t):
    chi_source = chi_magics_physical(energy_ev, bperp_t)
    energy_j = np.asarray(energy_ev, dtype=float) * ELEMENTARY_CHARGE
    scale = ALPHA * ELECTRON_MASS * C_LIGHT**2 / HBAR
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        k = kv(1.0 / 3.0, 2.0 / (3.0 * chi_source))
        return COEFFICIENT_ERBER * scale * (ELECTRON_MASS * C_LIGHT**2 / energy_j) * k**2


def build_grid(energy_ev=None, bperp_t=None, n_energy=25, n_b=25):
    if energy_ev is not None and bperp_t is not None:
        pairs = [(float(energy_ev), float(bperp_t), "requested")]
    else:
        energies = np.logspace(18, 21, n_energy)
        fields = np.logspace(-7, -4, n_b)
        pairs = [(float(e), float(b), "grid") for e in energies for b in fields]
    pairs.append((7.0e19, 2.115138e-5, "diagnostic"))
    return pd.DataFrame(pairs, columns=["E_eV", "Bperp_T", "case"])


def compare(output_dir, energy_ev=None, bperp_t=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    df = build_grid(energy_ev, bperp_t)
    chi_m = chi_magics_physical(df["E_eV"].to_numpy(), df["Bperp_T"].to_numpy())
    rate_m = magics_rate(df["E_eV"].to_numpy(), df["Bperp_T"].to_numpy())
    rate_e = erber1966_rate(df["E_eV"].to_numpy(), df["Bperp_T"].to_numpy())
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = rate_m / rate_e
        rel = np.abs(rate_m - rate_e) / np.maximum(np.abs(rate_e), np.finfo(float).tiny)

    out = pd.DataFrame(
        {
            "case": df["case"],
            "source": "Erber1966_Eq3.4_approx",
            "E_eV": df["E_eV"],
            "Bperp_T": df["Bperp_T"],
            "chi_M": chi_m,
            "chi_source": chi_m,
            "rate_MaGICS": rate_m,
            "rate_source": rate_e,
            "ratio": ratio,
            "relative_difference": rel,
        }
    )
    csv_path = output_dir / "erber_convention_comparison.csv"
    out.to_csv(csv_path, index=False)

    finite_ratio = out["ratio"].to_numpy()
    finite_ratio = finite_ratio[np.isfinite(finite_ratio)]
    diagnostic = out[out["case"] == "diagnostic"].iloc[-1]
    scale = ALPHA * ELECTRON_MASS * C_LIGHT**2 / HBAR
    c_m = RATE_PREFACTOR_MAGICS / scale
    ratio_prefactor = c_m / COEFFICIENT_ERBER

    summary = output_dir / "erber_convention_summary.md"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# Erber Convention Numerical Comparison\n\n")
        handle.write("Reference evaluated: Erber 1966 Eq. (3.4) using the analytic approximation Eq. (3.3d), with the photon parameter `x = chi_M`.\n\n")
        handle.write(f"- MaGICS coefficient `C_M`: `{c_m:.17e}`\n")
        handle.write(f"- Erber coefficient `C_E`: `{COEFFICIENT_ERBER:.17e}`\n")
        handle.write(f"- Prefactor ratio `C_M / C_E`: `{ratio_prefactor:.17e}`\n")
        if finite_ratio.size:
            handle.write(f"- Finite numerical ratio range: [`{finite_ratio.min():.17e}`, `{finite_ratio.max():.17e}`]\n")
        handle.write("\n## Diagnostic Case\n\n")
        for column in [
            "E_eV",
            "Bperp_T",
            "chi_M",
            "chi_source",
            "rate_MaGICS",
            "rate_source",
            "ratio",
            "relative_difference",
        ]:
            handle.write(f"- `{column}`: `{diagnostic[column]:.17e}`\n")
        handle.write("\nThe ratio is constant wherever both rates are finite because the Bessel argument and energy dependence are algebraically identical; only the dimensionless prefactor differs.\n")

    plt.figure()
    finite = np.isfinite(out["ratio"]) & (out["chi_M"] > 0)
    plt.semilogx(out.loc[finite, "chi_M"], out.loc[finite, "ratio"], ".")
    plt.xlabel("chi_M")
    plt.ylabel("rate_MaGICS / rate_Erber1966")
    plt.tight_layout()
    plt.savefig(output_dir / "erber_convention_ratio.png", dpi=160)
    plt.close()

    return {
        "rows": len(out),
        "comparison_csv": str(csv_path),
        "summary": str(summary),
        "prefactor_ratio": ratio_prefactor,
        "diagnostic_ratio": float(diagnostic["ratio"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("diagnostics"))
    parser.add_argument("--energy-eV", type=float, default=None)
    parser.add_argument("--Bperp-T", type=float, default=None)
    args = parser.parse_args()
    result = compare(args.output_dir, args.energy_eV, args.Bperp_T)
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
