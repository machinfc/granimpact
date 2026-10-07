"""Post-processing CLI: python3 -m granimpact_post <command> [...]"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .results import Case, load_case, load_many, sweep_table

SUBCOMMANDS = ("summary", "list", "sweep", "plot", "fit")


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
        print(f"\ntabla escrita en {out}")
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

    p = sub.add_parser("fit", help="fit D(E) to the analytical models")
    p.add_argument("inputs", nargs="+")
    p.set_defaults(func=cmd_fit)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
