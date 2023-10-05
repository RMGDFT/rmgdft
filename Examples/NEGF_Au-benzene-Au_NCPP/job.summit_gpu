#!/bin/bash
# Begin LSF directives
#BSUB -P MAT201
##BSUB -P MAT151
#BSUB -J N72
#BSUB -o tst.o%J
#BSUB -W 1:59
#BSUB -nnodes 1
#BSUB -alloc_flags smt4
#BSUB -alloc_flags gpumps
# End LSF directives and begin shell commands
module load gcc
module load boost
module load openblas/0.3.15-omp
module load cuda/11.0.3
module load fftw
module load hdf5
module load cmake/3.20.2
module load bzip2

#ldd ./rmg-gpu
export OMP_NUM_THREADS=1
export RMG_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_WAIT_POLICY=passive
export PAMI_IBV_ENABLE_DCT=1
export PAMI_ENABLE_STRIPING=0
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_IBV_DEVICE_NAME=mlx5_0:1
export PAMI_IBV_DEVICE_NAME_1=mlx5_3:1

export rmgon="/gpfs/alpine/proj-shared/mat151/rmg_exec/rmg-on-gpu-summit"
export rmgnegf="/gpfs/alpine/proj-shared/mat151/rmg_exec/rmg-negf-gpu-summit"
#jsrun -n144 -a1 -g1 -r6 -c7 --bind none --latency_priority cpu-memory --smpiargs "-gpu" ./rmg-gpu input.scalapack
cd lead1
jsrun -n6 -a1 -g1 -r6 -c1 --bind none --latency_priority cpu-memory --smpiargs "-gpu" $rmgon input
cd ../center
jsrun -n6 -a1 -g1 -r6 -c1 --bind none --latency_priority cpu-memory --smpiargs "-gpu" $rmgon input
cd ../3lead_lead1
jsrun -n6 -a1 -g1 -r6 -c1 --bind none --latency_priority cpu-memory --smpiargs "-gpu" $rmgnegf input
cd ../bias_0.0
jsrun -n6 -a1 -g1 -r6 -c1 --bind none --latency_priority cpu-memory --smpiargs "-gpu" $rmgnegf input
jsrun -n6 -a1 -g1 -r6 -c1 --bind none --latency_priority cpu-memory --smpiargs "-gpu" $rmgnegf input.110

