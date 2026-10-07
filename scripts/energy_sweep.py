#!/usr/bin/env python3
"""Energy sweep: D vs E, one case per speed.

Reuses `configs/energy_sweep.json` (the `sweep` block) and launches one case per value
changing only the projectile speed; at the end it fits the measured data to the
two models compared in this work (power law and logarithmic law).

    ./scripts/energy_sweep.py --dry-run    # see what would be launched
    ./scripts/energy_sweep.py              # launch everything (hours in serial)
    ./scripts/energy_sweep.py --only 2.0 4.0     # subset of speeds
    ./scripts/energy_sweep.py --no-vtk --jobs 1  # faster cases, one at a time
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "energy_sweep.json"
sys.path.insert(0, str(ROOT / "python"))


def find_binary() -> Path:
    for candidate in (
        ROOT / "build" / "omp" / "bin" / "granimpact_openmp",
        ROOT / "build" / "release" / "bin" / "granimpact_serial",
        ROOT / "build" / "cuda" / "bin" / "granimpact_cuda",
        ROOT / "build" / "debug" / "bin" / "granimpact_serial",
    ):
        if candidate.is_file():
            return candidate
    raise SystemExit("cannot find the executable; run ./scripts/install.sh")


def case_name(speed: float) -> str:
    return f"energy_sweep_v{speed:g}".replace(".", "p")


def run_case(binary: Path, speed: float, extra: list[str], verbose: bool) -> dict:
    name = case_name(speed)
    out_dir = ROOT / "data" / "outputs" / name
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(CONFIG), "--out", str(out_dir), "--name", name,
           "--speed", f"{speed:g}", *extra]
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start
    if verbose and result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        print(result.stderr or result.stdout, file=sys.stderr)
    return {"case": name, "v": speed, "code": result.returncode, "seconds": round(elapsed, 1),
            "output": str(out_dir)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", nargs="+", type=float, help="velocidades concretas [m/s]")
    parser.add_argument("--jobs", type=int, default=1, help="cases in parallel (default 1)")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="extra arguments for the simulator")
    parser.add_argument("--no-vtk", action="store_true", help="do not write VTK (faster)")
    parser.add_argument("--dry-run", action="store_true", help="only show the commands")
    parser.add_argument("--no-fit", action="store_true", help="do not fit the models at the end")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    speeds = config.get("sweep", {}).get("values", [])
    if not speeds:
        raise SystemExit(f"{CONFIG} has no sweep.values block")
    if args.only:
        speeds = [v for v in speeds if any(abs(v - o) < 1e-9 for o in args.only)]
        if not speeds:
            raise SystemExit("none of the requested speeds is in the config")

    binary = find_binary()
    extra = list(args.extra)
    if args.no_vtk and "--no-vtk" not in extra:
        extra.append("--no-vtk")

    print(f"Binario : {binary}")
    print(f"Config  : {CONFIG}")
    print(f"Cases   : {len(speeds)} -> {', '.join(f'{v:g} m/s' for v in speeds)}")
    print(f"Outputs : {ROOT / 'data' / 'outputs'}")

    if args.dry_run:
        for speed in speeds:
            out_dir = ROOT / "data" / "outputs" / case_name(speed)
            print("  " + " ".join([str(binary), str(CONFIG), "--out", str(out_dir),
                                   "--name", case_name(speed), "--speed", f"{speed:g}", *extra]))
        return 0

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        rows = list(pool.map(lambda s: run_case(binary, s, extra, args.verbose), speeds))

    print("\n=== Cases finished ===")
    for row in rows:
        status = "OK" if row["code"] == 0 else f"FAILED (code {row['code']})"
        print(f" {row['case']:<22} v={row['v']:>5g} m/s {status} {row['seconds']:>8.1f} s")

    failed = [r for r in rows if r["code"] != 0]
    if failed:
        return 1

    if not args.no_fit:
        try:
            from granimpact_post.models import compare_models
            from granimpact_post.results import load_case
        except ImportError as exc:  # pragma: no cover
            print(f"\n(could not import granimpact_post: {exc})")
            return 0
        cases = []
        for row in rows:
            csv_path = Path(row["output"]) / f"{row['case']}_crater_evolution.csv"
            if csv_path.is_file():
                cases.append(load_case(csv_path))
        if len(cases) < 2:
            print("\n(2+ cases with a CSV are needed to fit the models)")
            return 0
        print("\n=== D(E) fit with the measured data ===")
        result = compare_models(cases)
        for name in ("power", "logarithm"):
            print(f"  {name}: {result[name]}")
        if "best_by_R2" in result:
            print(f" best by R2: {result['best_by_R2']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
