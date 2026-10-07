#!/usr/bin/env bash
# Quick case: walks the full pipeline (bed -> settling -> impact -> analysis)
# with the `configs/default.json` case and leaves the results in data/outputs.
#
# ./scripts/run_quick.sh              # full case from the config
# ./scripts/run_quick.sh --no-vtk     # without VTK files (faster)
# ./scripts/run_quick.sh --t-end 0.02 # shorter impact
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN=""
BUILD_ROOT="${GRANIMPACT_BUILD_ROOT:-$ROOT/build}"   # see install.sh --build-dir
for c in "$BUILD_ROOT/release/bin/granimpact_serial" "$BUILD_ROOT/debug/bin/granimpact_serial"; do
    [[ -x "$c" ]] && BIN="$c" && break
done
[[ -n "$BIN" ]] || { echo "ERROR: build first with ./scripts/install.sh" >&2; exit 127; }
exec "$BIN" "$ROOT/configs/default.json" --out "$ROOT/data/outputs" "$@"
