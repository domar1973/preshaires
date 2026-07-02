#!/usr/bin/env python3
"""Analyze simple algebraic candidates for MaGICS pair prefactor."""

import argparse
import math
from fractions import Fraction
from pathlib import Path

import pandas as pd


C_LIGHT = 299792458.0
HBAR = 1.05457e-34
ELECTRON_MASS = 9.109534e-31
ALPHA = 1.0 / 137.0
MAGICS_PREFACTOR = 1.234e18
ERBER_COEFFICIENT = 0.16


def scale():
    return ALPHA * ELECTRON_MASS * C_LIGHT**2 / HBAR


def add_candidate(rows, name, value, formula, evidence="numeric candidate only"):
    target = MAGICS_PREFACTOR / scale()
    rows.append(
        {
            "name": name,
            "formula": formula,
            "value": value,
            "absolute_difference": abs(value - target),
            "relative_difference": abs(value - target) / abs(target),
            "prefactor_s^-1": value * scale(),
            "evidence": evidence,
        }
    )


def candidate_rows():
    rows = []
    target = MAGICS_PREFACTOR / scale()
    add_candidate(rows, "MaGICS", target, "1.234e18 / (alpha m_e c^2 / hbar)", "direct code literal")
    add_candidate(rows, "Erber Eq. 3.4", ERBER_COEFFICIENT, "0.16", "local primary Erber 1966")

    for n in range(1, 17):
        add_candidate(rows, f"sqrt(3)/{n}", math.sqrt(3) / n, f"sqrt(3)/{n}")
        add_candidate(rows, f"{n}*sqrt(3)/(8*pi)", n * math.sqrt(3) / (8 * math.pi), f"{n}*sqrt(3)/(8*pi)")
        add_candidate(rows, f"{n}*sqrt(3)/(16*pi)", n * math.sqrt(3) / (16 * math.pi), f"{n}*sqrt(3)/(16*pi)")
        add_candidate(rows, f"0.16*{n}/pi", ERBER_COEFFICIENT * n / math.pi, f"0.16*{n}/pi")
        add_candidate(rows, f"0.16*pi/{n}", ERBER_COEFFICIENT * math.pi / n, f"0.16*pi/{n}")
        add_candidate(rows, f"0.16*sqrt(3)*{n}/8", ERBER_COEFFICIENT * math.sqrt(3) * n / 8, f"0.16*sqrt(3)*{n}/8")
        add_candidate(rows, f"0.16*sqrt(3)*{n}/pi", ERBER_COEFFICIENT * math.sqrt(3) * n / math.pi, f"0.16*sqrt(3)*{n}/pi")

    for p in [-2, -1, 0, 1, 2]:
        for q in [-2, -1, 0, 1, 2]:
            for frac in [Fraction(a, b) for b in range(1, 17) for a in range(1, 33)]:
                value = float(frac) * (math.sqrt(3) ** p) * (math.pi ** q)
                add_candidate(
                    rows,
                    f"{frac}*sqrt(3)^{p}*pi^{q}",
                    value,
                    f"{frac} * sqrt(3)^{p} * pi^{q}",
                )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("diagnostics/pair_prefactor_candidates.csv"))
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    df = pd.DataFrame(candidate_rows())
    df = df.sort_values(["relative_difference", "formula"]).drop_duplicates("formula")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    target = MAGICS_PREFACTOR / scale()
    print(f"scale_s^-1: {scale():.17e}")
    print(f"target_C_M: {target:.17e}")
    print(f"target_over_Erber_0.16: {target / ERBER_COEFFICIENT:.17e}")
    print(f"output: {args.output}")
    print()
    print(df.head(args.top).to_string(index=False))


if __name__ == "__main__":
    main()
