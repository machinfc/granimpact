"""Post-processing tests (standard library, no dependencies)."""
from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from granimpact_post import load_case, sweep_table  # noqa: E402
from granimpact_post.mass import projectile_mass_kg  # noqa: E402
from granimpact_post.models import compare_models  # noqa: E402

CSV = """# GranImpact default | dt=5e-06 t_settle=0.25 t_end=0.08 phi=0.524 v=2.0 angle=90 gz=-9.81
t,particles,contacts,kinetic_energy,potential_energy,max_speed,d_exc,D,epsilon,h_rim,aspect_zd,v_in,ejecta_mass,classification
0.000000,100,0,0.200000,10.000,2.0,,,,,,,0,No data
0.001000,100,12,0.150000,10.000,1.5,0.005,0.040,0.9,0.002,0.125,0.0000010,0.0,Simple
0.002000,100,30,0.010000,9.900,0.3,0.008,0.050,0.8,0.003,0.160,0.0000020,0.0,Simple
# result: case=default launch_speed_m_s=2.0 impact_speed_m_s=1.31 settled_phi=0.549 D=0.05 launch_energy_J=0.2
"""

SWEEP = {
    "energy_sweep_v1": 0.039,
    "energy_sweep_v2": 0.156,
    "energy_sweep_v4": 0.625,
}


class ResultsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.default = self.dir / "default_crater_evolution.csv"
        self.default.write_text(CSV, encoding="utf-8")
        # Three sweep cases with D ~ E^0.25 (Uehara law) for the fit.
        for name, energy in SWEEP.items():
            d_mm = 100.0 * energy**0.25
            path = self.dir / f"{name}_crater_evolution.csv"
            path.write_text(
                "# GranImpact {n} | dt=2e-06 t_settle=0.25 t_end=0.1 phi=0.524 v=2 angle=90 gz=-9.81\n"
                "t,particles,contacts,kinetic_energy,potential_energy,max_speed,d_exc,D,epsilon,"
                "h_rim,aspect_zd,v_in,ejecta_mass,classification\n"
                "0.000000,50000,0,{e},1.0,2.0,,,,,,,0,No data\n"
                "0.050000,50000,900,0.5,1.0,1.0,0.01,{d},0.9,0.003,0.2,0.001,0.0,Simple\n".format(
                    n=name, e=energy, d=d_mm / 1000.0
                ),
                encoding="utf-8",
            )

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_reads_header_and_series(self) -> None:
        case = load_case(self.default)
        self.assertEqual(case.name, "default")
        self.assertEqual(case.meta["v"], "2.0")
        self.assertEqual(case.rows, 3)
        self.assertEqual(len(case.series("t")), 3)

    def test_crater_free_spans_do_not_enter_as_zero(self) -> None:
        case = load_case(self.default)
        self.assertIsNone(case.data["D"][0])
        self.assertEqual(len(case.finite("D")), 2)

    def test_summary_uses_last_row(self) -> None:
        summary = load_case(self.default).summary()
        self.assertAlmostEqual(summary["D_mm"], 50.0)
        self.assertAlmostEqual(summary["d_exc_mm"], 8.0)
        self.assertEqual(summary["class"], "Simple")

    def test_projectile_mass_comes_from_energy(self) -> None:
        # KE = 0.2 J with v = 2.0 m/s -> m = 2*KE/v^2 = 0.1 kg
        self.assertAlmostEqual(projectile_mass_kg(load_case(self.default)), 0.1)

    def test_final_block_carries_measured_metrics(self) -> None:
        case = load_case(self.default)
        self.assertEqual(case.meta["impact_speed_m_s"], "1.31")
        self.assertEqual(case.summary()["impact_speed_m_s"], "1.31")

    def test_impact_is_detected_at_first_row_with_contacts(self) -> None:
        self.assertEqual(load_case(self.default).impact_row(), 1)

    def test_sweep_recovers_expected_slope(self) -> None:
        cases = [load_case(p) for p in sorted(self.dir.glob("energy_sweep_*.csv"))]
        rows = sweep_table(cases)
        self.assertEqual(len(rows), 3)
        fit = compare_models(cases)["power"]
        # The synthetic data were generated with D ~ E^0.25: the fit must see it.
        self.assertAlmostEqual(fit["b"], 0.25, places=2)
        self.assertGreater(fit["R2"], 0.999)

    def test_compare_models_exposes_points_for_the_plot(self) -> None:
        # Regression: plot_sweep() reads result["points"]; a stale Spanish key
        # ("puntos") used to crash `granimpact-post sweep --plot` with KeyError.
        cases = [load_case(p) for p in sorted(self.dir.glob("energy_sweep_*.csv"))]
        result = compare_models(cases)
        self.assertIn("points", result)
        self.assertEqual(len(result["points"]["E_J"]), 3)
        self.assertEqual(len(result["points"]["D_mm"]), 3)
        self.assertEqual(result["power"]["b_expected_uehara"], 0.25)

    def test_launch_energy_prefers_the_csv_footer(self) -> None:
        # The writer records launch_energy_J in the footer; the table must use it
        # instead of re-deriving the mass, which fails for settled beds (the
        # particles already touch at t = 0, so no untouched row exists).
        case = load_case(self.default)
        self.assertEqual(case.meta["launch_energy_J"], "0.2")
        self.assertAlmostEqual(sweep_table([case])[0]["E_J"], 0.2)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_sweep_plot_needs_two_cases(self) -> None:
        from granimpact_post.plots import plot_sweep

        with self.assertRaises(SystemExit):
            plot_sweep([load_case(self.default)], self.dir / "one.png")

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_sweep_plot_writes_a_file(self) -> None:
        from granimpact_post.plots import plot_sweep

        cases = [load_case(p) for p in sorted(self.dir.glob("energy_sweep_*.csv"))]
        out = self.dir / "sweep.png"
        self.assertTrue(plot_sweep(cases, out).is_file())

    def test_csv_without_own_header(self) -> None:
        plain = self.dir / "plain.csv"
        plain.write_text("t,D\n0.0,0.01\n1.0,0.02\n", encoding="utf-8")
        case = load_case(plain)
        self.assertEqual(case.columns, ["t", "D"])
        self.assertEqual(case.rows, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
