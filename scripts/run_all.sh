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
BUILD_ROOT="${GRANIMPACT_BUILD_ROOT:-$ROOT/build}"   # see install.sh --build-dir
for c in "$BUILD_ROOT/release/bin/granimpact_serial" "$BUILD_ROOT/omp/bin/granimpact_openmp" \
         "$BUILD_ROOT/debug/bin/granimpact_serial"; do
    [[ -x "$c" ]] && BIN="$c" && break
done
[[ -n "$BIN" ]] || { echo "ERROR: build first with ./scripts/install.sh" >&2; exit 127; }

# Arguments: case names are anything that is not a flag. Flags that take a value
# must consume it here (--t-settle 0.02): otherwise the value would be taken for
# a case name. `--` stops flag parsing.
EXTRA=()
CASES=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-vtk) EXTRA+=(--no-vtk); shift ;;
        --dt|--t-end|--t-settle|--seed|--speed|--name|--backend|--out)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 needs a value" >&2; exit 2; }
            EXTRA+=("$1" "$2"); shift 2 ;;
        --)       shift; while [[ $# -gt 0 ]]; do CASES+=("$1"); shift; done ;;
        -*)       EXTRA+=("$1"); shift ;;
        *)        CASES+=("$1"); shift ;;
    esac
done
if [[ ${#CASES[@]} -eq 0 ]]; then
    CASES=(default loose_vertical compact_vertical loose_oblique45 lunar_impact energy_sweep)
fi

printf 'Binary: %s\n' "$BIN"
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
        printf '%-16s %s\n' "$case_name" "$(grep -E '^  (D \(diameter\)|d_exc)' "$summary" | tr '\n' ' ' | tr -s ' ')"
    fi
done
[[ ${#failed[@]} -eq 0 ]] && echo "All cases finished." || { echo "Failed: ${failed[*]}"; exit 1; }
