#!/usr/bin/env bash
# =============================================================================
# Local mirror of the GitHub CI (.github/workflows/ci.yml).
#
# Run it before pushing: it repeats, on your machine, exactly what GitHub will
# run. The clang part is the one that macOS exercises in CI: if you have clang
# installed, this script catches the same failures without waiting for a push.
#
#   ./scripts/ci_local.sh                 # everything (a few minutes)
#   ./scripts/ci_local.sh --quick         # skip the installer job
#   ./scripts/ci_local.sh --compilers "g++ clang++"
#
# Jobs mirrored:
#   build      gcc/clang x Debug/Release, warnings as errors, then ctest
#   python     post-processing unit tests + the D(E) fit check
#   installer  install.sh end to end into a temporary prefix + a real short run
#
# Build trees go to $TMPDIR/granimpact-ci, never inside the project.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/granimpact-ci"
CACHE="$WORK/deps-cache"          # FetchContent downloads once per machine
SKIP_INSTALLER=0
COMPILERS="g++ clang++"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --quick)       SKIP_INSTALLER=1; shift ;;
        --compilers)   COMPILERS="${2:?missing list}"; shift 2 ;;
        -h|--help)     sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

JOBS="${CMAKE_BUILD_PARALLEL_LEVEL:-$( (command -v nproc >/dev/null && nproc) || echo 4)}"
mkdir -p "$WORK" "$CACHE"

FAILURES=()
step() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
fail() { printf '    \033[31mFAILED:\033[0m %s\n' "$*"; FAILURES+=("$*"); }
ok()   { printf '    \033[32mOK\033[0m\n'; }

# ------------------------------------------------------------------- build ----
for preset in debug release; do
    if [[ "$preset" == debug ]]; then build_type=Debug; else build_type=Release; fi
    for cxx in $COMPILERS; do
        if ! command -v "$cxx" >/dev/null 2>&1; then
            step "build · $preset · $cxx (skipped: not installed)"
            continue
        fi
        step "build · $preset · $cxx"
        dir="$WORK/$preset-${cxx%%+*}"
        gen=()
        command -v ninja >/dev/null 2>&1 && gen=(-G Ninja)
        if cmake -S "$ROOT" -B "$dir" "${gen[@]}" \
                -DCMAKE_BUILD_TYPE="$build_type" \
                -DCMAKE_CXX_COMPILER="$cxx" \
                -DGRANIMPACT_WARNINGS_AS_ERRORS=ON \
                -DFETCHCONTENT_BASE_DIR="$CACHE" >"$dir-configure.log" 2>&1 \
           && cmake --build "$dir" -j "$JOBS" >"$dir-build.log" 2>&1 \
           && ctest --test-dir "$dir" --output-on-failure -j "$JOBS" >"$dir-ctest.log" 2>&1; then
            ok
        else
            grep -E "error|warning" "$dir-build.log" | head -10 || true
            grep -E "Failed|failed" "$dir-ctest.log" | head -10 || true
            fail "build · $preset · $cxx (logs: $dir-{configure,build,ctest}.log)"
        fi
    done
done

# ------------------------------------------------------------------ python ----
step "python · post-processing tests"
if command -v python3 >/dev/null 2>&1; then
    if PYTHONPATH="$ROOT/python" python3 -m unittest discover -s "$ROOT/python/tests" >"$WORK/python-tests.log" 2>&1; then
        ok
    else
        tail -15 "$WORK/python-tests.log"
        fail "python · unit tests"
    fi
    step "python · D(E) fit check"
    if PYTHONPATH="$ROOT/python" python3 -c "
from granimpact_post.models import fit_power_law
energies = [0.04, 0.16, 0.64]
diameters = [100.0 * e ** 0.25 for e in energies]
fit = fit_power_law(energies, diameters)
assert abs(fit['b'] - 0.25) < 1e-6, fit
assert abs(fit['a'] - 100.0) < 1e-6, fit
assert fit['R2'] > 0.999999, fit
print('fit OK: b =', fit['b'], '| R2 =', fit['R2'])
" ; then
        ok
    else
        fail "python · fit check"
    fi
else
    step "python (skipped: no python3)"
fi

# --------------------------------------------------------------- installer ----
if [[ $SKIP_INSTALLER -eq 1 ]]; then
    step "installer (skipped: --quick)"
else
    step "installer · install.sh end to end"
    PREFIX="$WORK/prefix"
    rm -rf "$PREFIX"
    if "$ROOT/scripts/install.sh" --prefix "$PREFIX" --jobs "$JOBS" >"$WORK/installer.log" 2>&1; then
        export PATH="$PREFIX/bin:$PATH"
        OUT="$WORK/out"
        rm -rf "$OUT"
        if granimpact --backends >/dev/null 2>&1 \
           && granimpact "$ROOT/configs/default.json" --out "$OUT" \
                     --t-settle 0.002 --t-end 0.0005 --no-vtk >"$WORK/installer-run.log" 2>&1 \
           && [[ -s "$OUT/default_summary.txt" && -s "$OUT/default_crater_evolution.csv" ]] \
           && PYTHONPATH="$ROOT/python" python3 -m granimpact_post summary "$OUT/default_crater_evolution.csv" >/dev/null 2>&1; then
            ok
        else
            tail -10 "$WORK/installer-run.log" || true
            fail "installer · short run / post-processing"
        fi
    else
        tail -15 "$WORK/installer.log"
        fail "installer · install.sh (log: $WORK/installer.log)"
    fi
fi

# ----------------------------------------------------------------- summary ----
if [[ ${#FAILURES[@]} -eq 0 ]]; then
    printf '\n\033[32mAll CI checks passed locally.\033[0m Safe to push.\n'
else
    printf '\n\033[31m%d check(s) failed:\033[0m\n' "${#FAILURES[@]}"
    printf '  - %s\n' "${FAILURES[@]}"
    exit 1
fi
