#!/usr/bin/env bash
# Runs the six project cases in series, one after another, saving each one
# in its own output folder. Meant to be left running overnight: it prints the time of
# each case and a final summary.
#
#   ./scripts/run_all.sh                 # all of them
#   ./scripts/run_all.sh default loose_vertical
# ./scripts/run_all.sh --no-vtk     # no VTK (much faster, no animation)
#
# COST WARNING: `loose_vertical`, `compact_vertical`, `loose_oblique45`,
# `lunar_impact` and every point of `energy_sweep` have ~50k particles and dt of
# 2e-6 s: hours per case in serial. That is exactly why stage 4 exists
# (OpenMP). See docs/roadmap.md.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN=""
for c in "$ROOT/build/release/bin/granimpact_serial" "$ROOT/build/omp/bin/granimpact_openmp" \
         "$ROOT/build/debug/bin/granimpact_serial"; do
    [[ -x "$c" ]] && BIN="$c" && break
done
[[ -n "$BIN" ]] || { echo "ERROR: build first with ./scripts/install.sh" >&2; exit 127; }

EXTRA=()
CASES=()
for arg in "$@"; do
    case "$arg" in
        --no-vtk) EXTRA+=(--no-vtk) ;;
        -*)       EXTRA+=("$arg") ;;
        *)        CASES+=("$arg") ;;
    esac
done
if [[ ${#CASES[@]} -eq 0 ]]; then
    CASES=(default loose_vertical compact_vertical loose_oblique45 lunar_impact energy_sweep)
fi

printf 'Binario: %s\n' "$BIN"
failed=()
for case_name in "${CASES[@]}"; do
    cfg="$ROOT/configs/$case_name.json"
    [[ -f "$cfg" ]] || { echo " ! $cfg does not exist"; failed+=("$case_name"); continue; }
    out="$ROOT/data/outputs/$case_name"
    mkdir -p "$out"
    printf '\n=== %s -> %s ===\n' "$case_name" "$out"
    start=$SECONDS
    if "$BIN" "$cfg" --out "$out" "${EXTRA[@]}"; then
        printf '    (%d s)\n' "$((SECONDS - start))"
    else
        echo " ! failed"
        failed+=("$case_name")
    fi
done

printf '\n=== Summary ===\n'
for case_name in "${CASES[@]}"; do
    summary="$ROOT/data/outputs/$case_name/${case_name}_summary.txt"
    if [[ -f "$summary" ]]; then
        printf '%-16s %s\n' "$case_name" "$(grep -E '^  (D \(diameter\)|d_exc) "$summary" | tr '\n' ' ' | tr -s ' ')"
    fi
done
[[ ${#failed[@]} -eq 0 ]] && echo "All cases finished." || { echo "Failed: ${failed[*]}"; exit 1; }
