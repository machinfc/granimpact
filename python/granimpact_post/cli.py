"""Post-processing CLI: python3 -m granimpact_post <command> [...]"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .results import (
    Case,
    export_craterslab_npz,
    load_case,
    load_many,
    resolve_sibling_csv,
    sweep_table,
)

SUBCOMMANDS = (
    "summary",
    "list",
    "sweep",
    "plot",
    "fit",
    "chains",
    "rim-theory",
    "formation",
    "visualize-all",
    "export-npz",
)


def _expand(patterns: list[str]) -> list[Path]:
    """Accepts paths and globs (expanded by the shell or by us)."""
    paths: list[Path] = []
    for pattern in patterns:
        path = Path(pattern)
        if path.is_file():
            paths.append(path)
        else:
            paths.extend(sorted(p for p in Path().glob(pattern) if p.is_file()))
    if not paths:
        raise SystemExit("no CSV files match those patterns")
    return paths


def _print_table(rows: list[dict], out=sys.stdout) -> None:
    if not rows:
        print("(no rows)", file=out)
        return
    columns = list(rows[0])
    widths = {c: max(len(c), *(len(_fmt(r[c])) for r in rows)) for c in columns}
    line = "  ".join(c.ljust(widths[c]) for c in columns)
    print(line, file=out)
    print("  ".join("-" * widths[c] for c in columns), file=out)
    for row in rows:
        print("  ".join(_fmt(row[c]).ljust(widths[c]) for c in columns), file=out)


def _fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def cmd_summary(args) -> int:
    for path in _expand(args.inputs):
        case = load_case(path)
        print(f"\n== {case.name}  ({path}) ==")
        print(case.as_text())
    return 0


def cmd_list(args) -> int:
    rows = []
    for path in _expand(args.inputs):
        case = load_case(path)
        rows.append(
            {
                "case": case.name,
                "rows": case.rows,
                "phi": case.meta.get("phi"),
                "v": case.meta.get("v"),
                "angle": case.meta.get("angle"),
                "path": str(path),
            }
        )
    _print_table(rows, sys.stdout)
    return 0


def cmd_sweep(args) -> int:
    cases = load_many(_expand(args.inputs))
    rows = sweep_table(cases)
    rows.sort(key=lambda r: (r["E_J"] is None, r["E_J"]))
    _print_table(rows, sys.stdout)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\ntable written to {out}")
    if args.plot:
        from .plots import plot_sweep

        path = plot_sweep(cases, args.plot)
        print(f"plot written to {path}")
    return 0


def cmd_fit(args) -> int:
    from .models import compare_models

    cases = load_many(_expand(args.inputs))
    result = compare_models(cases)
    for name, fit in result.items():
        print(f"{name}: {fit}")
    return 0


def cmd_plot(args) -> int:
    from .plots import plot_case

    for path in _expand(args.inputs):
        case = load_case(path)
        target = args.out or f"{case.name}.png"
        print(f"plot written to {plot_case(case, target)}")
    return 0


def cmd_chains(args) -> int:
    from .plots import plot_subsurface_chains

    for path in _expand(args.inputs):
        case = load_case(path)
        sub_csv = resolve_sibling_csv(case, "subsurface_chains")
        target = args.out or (path.parent / f"{case.name}_force_chains.png")
        print(f"force-chains plot written to {plot_subsurface_chains(case, sub_csv, target)}")
    return 0


def cmd_rim_theory(args) -> int:
    from .plots import plot_lateral_saturation

    for path in _expand(args.inputs):
        case = load_case(path)
        target = args.out or (path.parent / f"{case.name}_lateral_saturation.png")
        print(f"lateral-saturation plot written to {plot_lateral_saturation(case, target)}")
    return 0


def cmd_formation(args) -> int:
    from .plots import plot_crater_formation

    for path in _expand(args.inputs):
        case = load_case(path)
        rad_csv = resolve_sibling_csv(case, "radial_evolution")
        dmap_csv = resolve_sibling_csv(case, "depth_map")
        target = args.out or (path.parent / f"{case.name}_crater_formation.png")
        print(f"crater-formation plot written to {plot_crater_formation(case, rad_csv, dmap_csv, target)}")
    return 0


def cmd_export_npz(args) -> int:
    for path in _expand(args.inputs):
        dmap_csv = (
            path
            if path.name.endswith("_depth_map.csv")
            else resolve_sibling_csv(path, "depth_map")
        )
        out_path = export_craterslab_npz(
            dmap_csv,
            out_npz=args.out,
            use_total_surface=not args.bed_only,
        )
        print(f"craterslab .npz written to {out_path}")
    return 0


def cmd_visualize_all(args) -> int:
    from .plots import (
        plot_case,
        plot_crater_formation,
        plot_lateral_saturation,
        plot_subsurface_chains,
    )

    for path in _expand(args.inputs):
        case = load_case(path)
        out_dir = Path(args.out_dir) if args.out_dir else path.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"1) overview: {plot_case(case, out_dir / f'{case.name}_overview.png')}")
        print(f"2) rim/lateral saturation: {plot_lateral_saturation(case, out_dir / f'{case.name}_lateral_saturation.png')}")
        sub_csv = resolve_sibling_csv(case, "subsurface_chains")
        print(f"3) subsurface force chains: {plot_subsurface_chains(case, sub_csv, out_dir / f'{case.name}_force_chains.png')}")
        rad_csv = resolve_sibling_csv(case, "radial_evolution")
        dmap_csv = resolve_sibling_csv(case, "depth_map")
        print(f"4) crater formation: {plot_crater_formation(case, rad_csv, dmap_csv, out_dir / f'{case.name}_crater_formation.png')}")
        npz_path = export_craterslab_npz(dmap_csv, out_dir / f"{case.name}_craterslab.npz")
        print(f"5) craterslab .npz: {npz_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="granimpact_post",
        description="GranImpact results post-processing (CSV files from the simulator).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("summary", help="summary of one or more cases")
    p.add_argument("inputs", nargs="+", help="CSV o glob")
    p.set_defaults(func=cmd_summary)

    p = sub.add_parser("list", help="quick table of cases")
    p.add_argument("inputs", nargs="+")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("sweep", help="comparison table by energy (sweep)")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="write the table to a CSV")
    p.add_argument("--plot", help="also write a PNG plot")
    p.set_defaults(func=cmd_sweep)

    p = sub.add_parser("plot", help="plot a single case")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="output PNG")
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser("chains", help="Visualization 1: 2D subsurface force chains & Heckel saturation")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="output PNG")
    p.set_defaults(func=cmd_chains)

    p = sub.add_parser("rim-theory", help="Visualization 2: lateral opening & rim height saturation")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="output PNG")
    p.set_defaults(func=cmd_rim_theory)

    p = sub.add_parser("formation", help="Visualization 3: radial crater formation & 2D ToF depth map")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="output PNG")
    p.set_defaults(func=cmd_formation)

    p = sub.add_parser("export-npz", help="export 2D depth map to craterslab-compatible .npz")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--out", help="output .npz path")
    p.add_argument("--bed-only", action="store_true", help="export bed cavity without projectile grains")
    p.set_defaults(func=cmd_export_npz)

    p = sub.add_parser("visualize-all", help="generate all 3 scientific visualizations + craterslab .npz")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-d", "--out-dir", help="output directory for generated plots and .npz")
    p.set_defaults(func=cmd_visualize_all)

    p = sub.add_parser("fit", help="fit D(E) to the analytical models")
    p.add_argument("inputs", nargs="+")
    p.set_defaults(func=cmd_fit)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
