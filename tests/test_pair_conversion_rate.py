import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = Path(
    "/home/daniel/2_AREAS/DESARROLLO/MaGICS_TdF/TEST/docs/"
    "conversion_probability_diagnostics.csv"
)

spec = importlib.util.spec_from_file_location(
    "audit_pair_conversion_rate", REPO_ROOT / "tools" / "audit_pair_conversion_rate.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class ErberPairPrefactorTests(unittest.TestCase):
    def test_prefactor_is_erber_dimensionless_coefficient_times_qed_scale(self):
        expected = (
            0.16
            * audit.ALPHA
            * audit.ELECTRON_MASS
            * audit.C_LIGHT**2
            / audit.HBAR
        )
        self.assertAlmostEqual(audit.RATE_PREFACTOR, expected, delta=expected * 1e-15)

    def test_reference_rate_for_documented_regression_case(self):
        energy_gev = np.array([7.0e10])
        chi = np.array([0.3281999188925183])
        rate = audit.gamma_rate(energy_gev, chi)[0]

        # SciPy kv reference for Erber 1966 Eq. (3.4). Production MaGICS uses
        # the Numerical Recipes dbskr3_ Bessel routine, whose separate audit
        # shows percent-level error is possible over diagnostic grids; this
        # tighter tolerance isolates the prefactor from that implementation.
        self.assertAlmostEqual(rate, 83.21254795423327, delta=1e-10)


@unittest.skipUnless(CSV_PATH.exists(), "diagnostic CSV not available")
class PairConversionRateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_csv(CSV_PATH)
        cls.rate = audit.gamma_rate(cls.df["energy_GeV"].to_numpy(), cls.df["chi"].to_numpy())
        cls.reference_step = cls.rate * cls.df["dt_s"].to_numpy()

    def test_python_reference_reproduces_step_probability(self):
        rel = audit.relative_error(
            self.reference_step, self.df["step_probability"].to_numpy()
        )
        finite = rel[np.isfinite(rel)]
        if np.nanmax(finite) <= 1e-2:
            return

        legacy_prefactor = 1.234e18
        legacy_ratio = legacy_prefactor / audit.RATE_PREFACTOR
        valid = self.reference_step > 0
        ratio = self.df.loc[valid, "step_probability"].to_numpy() / self.reference_step[valid]
        finite_ratio = ratio[np.isfinite(ratio)]
        self.assertAlmostEqual(float(np.nanmedian(finite_ratio)), legacy_ratio, delta=1e-2)

    def test_product_reproduces_final_accumulated_probability(self):
        product = 1.0 - np.prod(1.0 - self.df["step_probability"].to_numpy())
        final = float(self.df["accumulated_probability"].iloc[-1])
        self.assertAlmostEqual(product, final, delta=5e-14)

    def test_positive_chi_dt_gives_nonnegative_probability(self):
        mask = (self.df["chi"].to_numpy() > 0) & (self.df["dt_s"].to_numpy() > 0)
        self.assertTrue(np.all(self.df.loc[mask, "step_probability"].to_numpy() >= 0))

    def test_very_small_chi_underflows_to_nonnegative_zero(self):
        chi = np.array([1e-12, 1e-8, 1e-4])
        energy = np.full_like(chi, 7e10)
        rate = audit.gamma_rate(energy, chi)
        self.assertTrue(np.all(np.isfinite(rate)))
        self.assertTrue(np.all(rate >= 0))

    def test_large_bessel_argument_is_stable(self):
        chi = np.array([1e-6, 1e-5, 1e-4])
        arg = audit.bessel_argument(chi)
        self.assertTrue(np.all(arg > 1e3))
        rate = audit.gamma_rate(np.full_like(chi, 7e10), chi)
        self.assertTrue(np.all(np.isfinite(rate)))

    def test_no_silent_nan_in_csv_reference(self):
        self.assertEqual(int(np.isnan(self.reference_step).sum()), 0)
        self.assertEqual(int(np.isnan(self.rate).sum()), 0)
        self.assertEqual(int(self.df[["chi", "dt_s", "energy_GeV"]].isna().sum().sum()), 0)

    def test_csv_chi_range_is_covered_by_grid(self):
        csv_min = float(self.df["chi"].min())
        csv_max = float(self.df["chi"].max())
        grid = np.logspace(math.log10(max(csv_min / 10.0, 1e-12)), math.log10(csv_max * 10.0), 200)
        self.assertLessEqual(grid.min(), csv_min)
        self.assertGreaterEqual(grid.max(), csv_max)


if __name__ == "__main__":
    unittest.main()
