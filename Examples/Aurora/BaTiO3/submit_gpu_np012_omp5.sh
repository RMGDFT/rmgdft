#!/bin/bash
#PBS -A HeavyMetalQMC
#PBS -N rmg_gpu_np012_omp5
#PBS -l select=1
#PBS -l place=scatter
#PBS -l filesystems=flare
#PBS -q debug
#PBS -l walltime=00:10:00

set -u -o pipefail

echo "============================================"
echo "RMG-DFT PBE Benchmark"
echo "Platform : GPU"
echo "NP       : 12"
echo "OMP      : 5"
echo "kpdist   : 12"
echo "ProcGrid : 1x1x1"
echo "============================================"

cd /flare/HeavyMetalQMC/abenali/BaTiO3/bench/bench_inputs/run_gpu_np012_omp5 || exit 1
echo "Jobid: $PBS_JOBID"
echo "Node:"
cat $PBS_NODEFILE

# --- Modules ---
for m in oneapi/release oneapi/eng-compiler; do
  module is-loaded "$m" && module unload "$m"
done
module load oneapi/release/2025.2.0 mpich/opt/develop-git.6037a7a
module load boost hdf5 netlib-scalapack

# --- GPU Environment ---
export MPIR_CVAR_ENABLE_GPU=0
export FI_CXI_DEFAULT_CQ_SIZE=131072
export LIBOMP_USE_HIDDEN_HELPER_TASK=0
export ZES_ENABLE_SYSMAN=1
export ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE
export NEOReadDebugKeys=1
export SplitBcsCopy=0
export OMP_PLACES=cores
export OMP_NUM_THREADS=5

echo "Starting run at $(date)"
start=$(date +%s)

mpiexec --hostfile $PBS_NODEFILE \
        -np 12 -ppn 12 \
        -d 5 --cpu-bind depth \
        /soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh \
        /home/abenali/Software/rmgdft/build_aurora_oneapi2025.2.0_gpu/rmg-gpu BaTiO3_PBE_gpu_np012_omp5.in \
        > rmg_run.out 2>&1

EXIT_CODE=$?
end=$(date +%s)
ELAPSED=$((end - start))

echo "Exit code : $EXIT_CODE"
echo "Wall time : ${ELAPSED}s"

# --- Extract timing ---
echo ""
echo "--- SCF Timing ---"
grep -E "quench|elapsed|SCF" rmg_run.out | tail -10

# --- Append to global summary ---
SUMMARY=/flare/HeavyMetalQMC/abenali/BaTiO3/bench/bench_inputs/bench_summary_global.txt
echo "GPU |  12 | 5 |  12 | 1x1x1 | $ELAPSED | $EXIT_CODE" >> $SUMMARY
echo "Done at $(date)"
