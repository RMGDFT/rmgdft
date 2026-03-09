#!/bin/bash
# =============================================================================
# install_rmgdft_aurora.sh
# Build and install RMG-DFT on Aurora (Intel Data Center GPU Max / PVC)
# with oneAPI 2025.2.0 and SYCL GPU support.
#
# Usage:
#   bash install_rmgdft_aurora.sh [--gpu | --cpu | --both] [--clean]
#
# Options:
#   --gpu    Build GPU-enabled binary only (default)
#   --cpu    Build CPU-only binary only
#   --both   Build both GPU and CPU binaries
#   --clean  Remove existing build directories before building
#
# Author:  A. Benali (Argonne National Laboratory)
# Tested:  Aurora, oneAPI 2025.2.0, Intel Data Center GPU Max 1550
# =============================================================================

set -e
set -o pipefail

# ---------------------------------------------------------------------------
# COLOR OUTPUT
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; exit 1; }
header()  { echo -e "\n${BOLD}${CYAN}========== $* ==========${RESET}\n"; }

# ---------------------------------------------------------------------------
# DEFAULTS
# ---------------------------------------------------------------------------
BUILD_GPU=true
BUILD_CPU=false
CLEAN=false

# ---------------------------------------------------------------------------
# PARSE ARGUMENTS
# ---------------------------------------------------------------------------
for arg in "$@"; do
    case "$arg" in
        --gpu)   BUILD_GPU=true;  BUILD_CPU=false ;;
        --cpu)   BUILD_GPU=false; BUILD_CPU=true  ;;
        --both)  BUILD_GPU=true;  BUILD_CPU=true  ;;
        --clean) CLEAN=true ;;
        --help)
            sed -n '3,16p' "$0"
            exit 0
            ;;
        *) warn "Unknown argument: $arg (ignored)" ;;
    esac
done

# ---------------------------------------------------------------------------
# PATHS  — edit these if your layout differs
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${SCRIPT_DIR}"                           # repo root
INSTALL_PREFIX="${HOME}/Software/rmgdft_install"    # optional install prefix

BUILD_GPU_DIR="${SOURCE_DIR}/build_aurora_oneapi2025.2.0_gpu"
BUILD_CPU_DIR="${SOURCE_DIR}/build_aurora_oneapi2025.2.0_cpu"

NPROC=$(nproc 2>/dev/null || echo 16)   # parallel make jobs

# ---------------------------------------------------------------------------
# VERIFY WE ARE ON AN AURORA LOGIN/COMPUTE NODE
# ---------------------------------------------------------------------------
header "Environment Check"

if [[ ! -d /opt/aurora ]]; then
    warn "  /opt/aurora not found — are you on an Aurora node?"
fi

# ---------------------------------------------------------------------------
# MODULES
# ---------------------------------------------------------------------------
header "Loading Modules"

# Unload conflicting compiler modules first
for mod in oneapi/release oneapi/eng-compiler; do
    if module is-loaded "$mod" 2>/dev/null; then
        info "Unloading $mod"
        module unload "$mod"
    fi
done

module load oneapi/release/2025.2.0       || error "Failed to load oneapi/release/2025.2.0"
module load mpich/opt/develop-git.6037a7a || error "Failed to load mpich"
module load boost                          || error "Failed to load boost"
module load hdf5                           || error "Failed to load hdf5"
module load netlib-scalapack               || error "Failed to load netlib-scalapack"

success "Modules loaded:"
module list 2>&1 | grep -E "oneapi|mpich|boost|hdf5|scalapack" | sed 's/^/  /'

# ---------------------------------------------------------------------------
# COMPILER / ENVIRONMENT SANITY CHECK
# ---------------------------------------------------------------------------
header "Compiler Check"

CXX_COMPILER=$(which icpx 2>/dev/null || which dpcpp 2>/dev/null || echo "")
C_COMPILER=$(which icx 2>/dev/null || echo "")
FC_COMPILER=$(which ifx 2>/dev/null || echo "")

[[ -z "$CXX_COMPILER" ]] && error "icpx/dpcpp not found — oneAPI not loaded correctly"

info "C   compiler : $C_COMPILER"
info "C++ compiler : $CXX_COMPILER"
info "Fortran      : $FC_COMPILER"
info "MPI wrapper  : $(which mpicc 2>/dev/null || echo 'not found')"

# ---------------------------------------------------------------------------
# CMAKE FLAGS SHARED BY BOTH BUILDS
# ---------------------------------------------------------------------------
COMMON_CMAKE_FLAGS=(
    -DCMAKE_C_COMPILER=mpicc
    -DCMAKE_CXX_COMPILER=mpicxx
    -DCMAKE_Fortran_COMPILER=mpifort
    -DCMAKE_BUILD_TYPE=Release
    -DENABLE_PHDF5=ON             # parallel HDF5 output (flag name from CMakeLists)
    -DUSE_INTERNAL_SCALAPACK=OFF  # use the module-loaded netlib-scalapack, not bundled one
    -DRMG_LIBXC_ENABLED=OFF       # libXC causes build failures on Aurora with oneAPI 2025; default is OFF but explicit here
)

# ---------------------------------------------------------------------------
# GPU BUILD
# ---------------------------------------------------------------------------
build_gpu() {
    header "Building GPU Binary (SYCL / Intel PVC)"

    if $CLEAN && [[ -d "$BUILD_GPU_DIR" ]]; then
        info "Removing existing GPU build directory..."
        rm -rf "$BUILD_GPU_DIR"
    fi

    mkdir -p "$BUILD_GPU_DIR"
    cd "$BUILD_GPU_DIR"

    info "Running CMake for GPU build..."
    cmake "${SOURCE_DIR}" \
        "${COMMON_CMAKE_FLAGS[@]}" \
        -DRMG_SYCL_ENABLED=ON        \
        2>&1 | tee cmake_gpu.log

    info "Building GPU binary with ${NPROC} parallel jobs..."
    make -j"${NPROC}" 2>&1 | tee make_gpu.log

    # Verify the binary was produced
    if [[ -x "${BUILD_GPU_DIR}/rmg-gpu" ]]; then
        success "GPU binary built: ${BUILD_GPU_DIR}/rmg-gpu"
    elif [[ -x "${BUILD_GPU_DIR}/bin/rmg-gpu" ]]; then
        success "GPU binary built: ${BUILD_GPU_DIR}/bin/rmg-gpu"
    else
        error "GPU binary not found after build — check ${BUILD_GPU_DIR}/make_gpu.log"
    fi

    cd "${SOURCE_DIR}"
}

# ---------------------------------------------------------------------------
# CPU BUILD
# ---------------------------------------------------------------------------
build_cpu() {
    header "Building CPU Binary"

    if $CLEAN && [[ -d "$BUILD_CPU_DIR" ]]; then
        info "Removing existing CPU build directory..."
        rm -rf "$BUILD_CPU_DIR"
    fi

    mkdir -p "$BUILD_CPU_DIR"
    cd "$BUILD_CPU_DIR"

    info "Running CMake for CPU build..."
    cmake "${SOURCE_DIR}" \
        "${COMMON_CMAKE_FLAGS[@]}" \
        -DRMG_SYCL_ENABLED=OFF       \
        -DCMAKE_CXX_FLAGS="-O3 -xCORE-AVX512 -fp-model=precise" \
        2>&1 | tee cmake_cpu.log

    info "Building CPU binary with ${NPROC} parallel jobs..."
    make -j"${NPROC}" 2>&1 | tee make_cpu.log

    if [[ -x "${BUILD_CPU_DIR}/rmg-cpu" ]] || [[ -x "${BUILD_CPU_DIR}/bin/rmg-cpu" ]]; then
        success "CPU binary built successfully"
    else
        error "CPU binary not found after build — check ${BUILD_CPU_DIR}/make_cpu.log"
    fi

    cd "${SOURCE_DIR}"
}

# ---------------------------------------------------------------------------
# RUN BUILDS
# ---------------------------------------------------------------------------
$BUILD_GPU && build_gpu
$BUILD_CPU && build_cpu

# ---------------------------------------------------------------------------
# QUICK SMOKE TEST
# ---------------------------------------------------------------------------
header "Smoke Test"

if $BUILD_GPU; then
    EXE=$(ls "${BUILD_GPU_DIR}/rmg-gpu" "${BUILD_GPU_DIR}/bin/rmg-gpu" 2>/dev/null | head -1)
    if [[ -x "$EXE" ]]; then
        info "Checking shared library dependencies..."
        ldd "$EXE" 2>&1 | grep -i "not found" && \
            warn "Some libraries not found — check module environment at runtime" || \
            success "All shared libraries resolved"
    fi
fi

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
header "Build Summary"

echo -e "  Source directory  : ${SOURCE_DIR}"
echo -e "  oneAPI version    : 2025.2.0"
echo -e "  Target GPU        : Intel Data Center GPU Max 1550 (PVC)"
echo ""

if $BUILD_GPU; then
    EXE=$(ls "${BUILD_GPU_DIR}/rmg-gpu" "${BUILD_GPU_DIR}/bin/rmg-gpu" 2>/dev/null | head -1)
    echo -e "  GPU binary        : ${EXE:-NOT FOUND}"
fi
if $BUILD_CPU; then
    EXE=$(ls "${BUILD_CPU_DIR}/rmg-cpu" "${BUILD_CPU_DIR}/bin/rmg-cpu" 2>/dev/null | head -1)
    echo -e "  CPU binary        : ${EXE:-NOT FOUND}"
fi

echo ""
echo -e "  Recommended run settings (1 node, 12 GPU tiles):"
echo -e "    export OMP_NUM_THREADS=5"
echo -e "    export OMP_PLACES=cores"
echo -e "    NRANKS=12, --cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"
echo -e "    /soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh <exe> <input>"
echo ""

success "Done."

