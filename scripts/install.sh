#!/usr/bin/env bash
# =============================================================================
# GranImpact · full local installation
# -----------------------------------------------------------------------------
# Builds everything that exists in the project today: compiles the library and
# the executable, runs the tests, installs binaries/configs/docs/python, creates
# the data directories and (optionally) runs a quick case.
#
# It is idempotent: it can be rerun after every code change.
#
# Usage:
#   ./scripts/install.sh                     # Release build in ./build/release + tests
#   ./scripts/install.sh --prefix ~/.local   # install prefix
#   ./scripts/install.sh --quick             # also runs the quick case
#   ./scripts/install.sh --debug             # Debug build with verbose tests
#   ./scripts/install.sh --omp               # enables the OpenMP backend (stage 4)
#   ./scripts/install.sh --no-tests --no-install
#   ./scripts/install.sh --no-python         # skip installing the Python package
#
# Environment variables honored: CC, CXX, CMAKE_BUILD_PARALLEL_LEVEL.
# =============================================================================
set -euo pipefail

# ------------------------------------------------------------------- output ---
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    B=$'\033[1m'; G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; N=$'\033[0m'
else
    B=''; G=''; Y=''; R=''; N=''
fi
step() { printf '%s==>%s %s\n' "$B" "$N" "$*"; }
ok()   { printf '    %s✓%s %s\n' "$G" "$N" "$*"; }
warn() { printf '    %s!%s %s\n' "$Y" "$N" "$*"; }
die()  { printf '\n%sERROR:%s %s\n\n' "$R" "$N" "$*" >&2; exit 1; }

# ----------------------------------------------------------------- options ----
PREFIX="${PREFIX:-$HOME/.local}"
BUILD_TYPE="Release"
PRESET="release"
JOBS="${CMAKE_BUILD_PARALLEL_LEVEL:-}"
RUN_TESTS=1
DO_INSTALL=1
DO_QUICK=0
ENABLE_OMP=0
ENABLE_CUDA=0
PYTHON_PKG=1
BUILD_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)     PREFIX="${2:?missing directory}"; shift 2 ;;
        --build-dir)  BUILD_DIR="${2:?missing directory}"; shift 2 ;;
        --jobs|-j)    JOBS="${2:?missing number}"; shift 2 ;;
        --debug)      BUILD_TYPE="Debug"; PRESET="debug"; shift ;;
        --omp)        ENABLE_OMP=1; PRESET="omp"; shift ;;
        --cuda)       ENABLE_CUDA=1; PRESET="cuda"; shift ;;
        --quick)      DO_QUICK=1; shift ;;
        --no-tests)   RUN_TESTS=0; shift ;;
        --no-install) DO_INSTALL=0; shift ;;
        --no-python)  PYTHON_PKG=0; shift ;;
        -h|--help)
            sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) die "unknown option: $1 (use --help)" ;;
    esac
done

# --------------------------------------------------------------- proyecto -----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
[[ -f "$ROOT/CMakeLists.txt" ]] || die "cannot find CMakeLists.txt in $ROOT"
BUILD_DIR="${BUILD_DIR:-$ROOT/build/$PRESET}"
LOG_DIR="$ROOT/build/logs"
JOBS="${JOBS:-$( (command -v nproc >/dev/null && nproc) || (command -v sysctl >/dev/null && sysctl -n hw.ncpu) || echo 4)}"

step "GranImpact · local installation"
printf '    proyecto ....... %s\n' "$ROOT"
printf '    build .......... %s (%s)\n' "$BUILD_DIR" "$BUILD_TYPE"
printf '    prefix ......... %s\n' "$PREFIX"
printf '    jobs ........... %s\n' "$JOBS"

# -------------------------------------------------------------- requisitos ----
step "Checking requirements"
need() { command -v "$1" >/dev/null 2>&1 || die "$2"; }
need cmake "missing cmake (>= 3.20). Linux: 'sudo apt install cmake'; macOS: 'brew install cmake'"
CMAKE_VERSION="$(cmake --version | head -1 | awk '{print $3}')"
ok "cmake $CMAKE_VERSION"

if [[ -n "${CXX:-}" ]]; then CXX_BIN="$CXX"
elif command -v g++ >/dev/null 2>&1; then CXX_BIN=g++
elif command -v clang++ >/dev/null 2>&1; then CXX_BIN=clang++
else die "no C++ compiler (g++ or clang++). Linux: 'sudo apt install build-essential'"; fi
ok "compiler: $CXX_BIN $("$CXX_BIN" -dumpversion 2>/dev/null || echo '?')"

command -v ninja >/dev/null 2>&1 && ok "ninja (faster builds)" || warn "no ninja: falling back to Makefiles (slower)"

# FetchContent downloads nlohmann/json and GoogleTest the first time.
FETCH_OK=1
if [[ -d "$BUILD_DIR/_deps" ]]; then
    ok "dependencias ya descargadas (build existente)"
    FETCH_OK=0
else
    for tool in git curl; do
        command -v "$tool" >/dev/null 2>&1 || { warn "missing $tool (needed by FetchContent)"; FETCH_OK=0; }
    done
    if command -v git >/dev/null 2>&1; then
        if git ls-remote --exit-code https://github.com/nlohmann/json.git HEAD >/dev/null 2>&1; then
            ok "GitHub reachable (FetchContent can download)"
        else
            warn "GitHub is not responding: if the dependencies are not on the system, configuration will fail"
            warn "  alternatives: -DGRANIMPACT_USE_SYSTEM_DEPS=ON (with nlohmann-json3-dev and libgtest-dev)"
        fi
    fi
fi

# ------------------------------------------------------------------ build -----
mkdir -p "$LOG_DIR" "$BUILD_DIR"

step "Configurando"
CMAKE_ARGS=(
    -S "$ROOT" -B "$BUILD_DIR"
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    -DCMAKE_INSTALL_PREFIX="$PREFIX"
    -DGRANIMPACT_BUILD_TESTS="$([[ $RUN_TESTS -eq 1 ]] && echo ON || echo OFF)"
    -DGRANIMPACT_ENABLE_OPENMP="$([[ $ENABLE_OMP -eq 1 ]] && echo ON || echo OFF)"
    -DGRANIMPACT_ENABLE_CUDA="$([[ $ENABLE_CUDA -eq 1 ]] && echo ON || echo OFF)"
)
if command -v ninja >/dev/null 2>&1; then CMAKE_ARGS+=(-G Ninja); fi
if ! cmake "${CMAKE_ARGS[@]}" >"$LOG_DIR/configure-$PRESET.log" 2>&1; then
    tail -30 "$LOG_DIR/configure-$PRESET.log" >&2
    die "configuration failed (full log in $LOG_DIR/configure-$PRESET.log)"
fi
ok "configured"

step "Building ($JOBS jobs)"
if ! cmake --build "$BUILD_DIR" -j "$JOBS" >"$LOG_DIR/build-$PRESET.log" 2>&1; then
    grep -E "error|Error" "$LOG_DIR/build-$PRESET.log" | head -20 >&2 || tail -20 "$LOG_DIR/build-$PRESET.log" >&2
    die "build failed (full log in $LOG_DIR/build-$PRESET.log)"
fi
ok "binaries in $BUILD_DIR/bin"

# ------------------------------------------------------------------ tests -----
if [[ $RUN_TESTS -eq 1 ]]; then
    step "Tests (ctest)"
    if (cd "$BUILD_DIR" && ctest --output-on-failure -j "$JOBS" >"$LOG_DIR/ctest-$PRESET.log" 2>&1); then
        ok "$(grep -E 'tests passed' "$LOG_DIR/ctest-$PRESET.log" | tail -1)"
    else
        grep -E "Failed|Failure" -A3 "$LOG_DIR/ctest-$PRESET.log" | head -30 >&2
        die "there are failing tests (full log in $LOG_DIR/ctest-$PRESET.log)"
    fi
else
    warn "tests skipped (--no-tests)"
fi

# --------------------------------------------------------------- instalar -----
if [[ $DO_INSTALL -eq 1 ]]; then
    step "Installing into $PREFIX"
    cmake --install "$BUILD_DIR" >/dev/null
    ok "bin/granimpact_serial, bin/granimpact (launcher), share/granimpact/{configs,docs,python}"
else
    warn "installation skipped (--no-install)"
fi

# ------------------------------------------------------------------ python ----
if [[ $PYTHON_PKG -eq 1 ]]; then
    if command -v python3 >/dev/null 2>&1; then
        if python3 -c "import matplotlib" >/dev/null 2>&1; then
            ok "python3 $(python3 -c 'import platform;print(platform.python_version())') with matplotlib (plots available)"
        else
            warn "python3 without matplotlib: post-processing works, plots do not ('pip install matplotlib')"
        fi
        if command -v pip3 >/dev/null 2>&1 && [[ $DO_INSTALL -eq 1 ]]; then
            if pip3 install --quiet --user "$ROOT/python" >/dev/null 2>&1; then
                ok "post-processing package installed (python3 -m granimpact_post --help)"
            else
                warn "could not install the package with pip; use PYTHONPATH=$ROOT/python"
            fi
        fi
    else
        warn "no python3: post-processing is disabled"
    fi
fi

# -------------------------------------------------------------------- data ----
mkdir -p "$ROOT/data/outputs"
ok "outputs directory: $ROOT/data/outputs"

# -------------------------------------------------------------------- quick ---
BIN="$BUILD_DIR/bin/granimpact_serial"
if [[ $DO_QUICK -eq 1 ]]; then
    step "Quick case (validates the full pipeline)"
    "$BIN" "$ROOT/configs/default.json" --out "$ROOT/data/outputs"
    ok "outputs in $ROOT/data/outputs"
fi

# --------------------------------------------------------------------- end ----
cat <<EOF

${B}Done.${N} GranImpact is installed.

  Executable ....... $BIN
  (also at) ......... $PREFIX/bin/granimpact_serial
  Launcher ......... $PREFIX/bin/granimpact   (picks the best available backend)
  Configs .......... $PREFIX/share/granimpact/configs
  Documentation .... $PREFIX/share/granimpact/docs
  Build logs ....... $LOG_DIR

Next steps
  $BIN $ROOT/configs/default.json --out $ROOT/data/outputs      # quick case
  $BIN $ROOT/configs/loose_vertical.json --no-vtk                # ~50k grains
  python3 -m granimpact_post summary $ROOT/data/outputs/*.csv    # post-processing

If the prefix is not in PATH:
  export PATH="$PREFIX/bin:\$PATH"
  export PYTHONPATH="$ROOT/python:\${PYTHONPATH:-}"
EOF
