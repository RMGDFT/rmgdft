#!/bin/bash
# Begin LSF directives
##BSUB -P CHP107
#BSUB -P MAT151
#BSUB -J N72
#BSUB -o tst.o%J
#BSUB -W 0:59
#BSUB -nnodes 36
#BSUB -alloc_flags smt4
# End LSF directives and begin shell commands
module load gcc/7.4.0
module load boost
module load openblas/0.3.9-omp
module load cuda
module load fftw
module load hdf5

cd /gpfs/alpine/proj-shared/mat151/elbriggs/Water1000/

#ldd ./rmg-gpu
export OMP_NUM_THREADS=6
export RMG_NUM_THREADS=6
export OMP_DYNAMIC=FALSE
export OMP_WAIT_POLICY=passive
export PAMI_IBV_ENABLE_DCT=1
export PAMI_ENABLE_STRIPING=0
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_IBV_DEVICE_NAME=mlx5_0:1
export PAMI_IBV_DEVICE_NAME_1=mlx5_3:1

jsrun -n216 -a1 -g1 -r6 -c7 --bind none --latency_priority cpu-memory --smpiargs "-gpu" ./rmg-gpu input
