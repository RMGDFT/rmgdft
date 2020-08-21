#!/bin/bash
# Begin LSF directives
#BSUB -P CHP107
##BSUB -P MAT151
#BSUB -J N96
#BSUB -o tst.o%J
#BSUB -W 0:30
#BSUB -nnodes 96
#BSUB -alloc_flags smt4
# End LSF directives and begin shell commands
module load gcc/6.4.0
module load boost
module load essl
module load cuda
module load fftw
module load hdf5

cd $PROJWORK/chp107/elbriggs/NiO512/
cp ~/bin/rmg-gpu $PROJWORK/chp107/elbriggs/NiO512/

export OMP_NUM_THREADS=6
export RMG_NUM_THREADS=6
export OMP_DYNAMIC=TRUE
export OMP_WAIT_POLICY=passive
export PAMI_IBV_ENABLE_DCT=1
export PAMI_ENABLE_STRIPING=0
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_IBV_DEVICE_NAME=mlx5_0:1
export PAMI_IBV_DEVICE_NAME_1=mlx5_3:1

jsrun -n576 -a1 -g1 -r6 -c7 --bind none --latency_priority cpu-memory ./rmg-gpu input
