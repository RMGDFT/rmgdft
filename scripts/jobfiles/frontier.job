#!/bin/bash
# Change to your account
#SBATCH -A MAT151
# Job name
#SBATCH -J NiO64
#SBATCH -o %x-%j.out
# Time
#SBATCH -t 00:30:00
# queue
#SBATCH -p batch
#
# Number of frontier nodes to use
#SBATCH -N 4
#
# OMP num threads. Frontier reserves 8 of 64 cores on a node
# for the system. There are 8 logical GPUs per node so we use
# 8 MPI tasks/node with 7 OMP threads per node
export OMP_NUM_THREADS=7
#
# RMG threads. Max of 7 same as for OMP_NUM_THREADS but in some
# cases running with fewer may yield better performance because
# of cache effects.
export RMG_NUM_THREADS=5
export MPICH_OFI_NIC_POLICY=NUMA
# Don't change this
export MPICH_GPU_SUPPORT_ENABLED=0
#
# Load modules
module load PrgEnv-gnu/8.3.3
module load bzip2
module load boost/1.79.0-cxx17
module load cray-fftw
module load cray-hdf5-parallel
module load craype-accel-amd-gfx90a
module load rocm/5.4.3
# ntasks should be 8*nodes where nodes is set with -N above
srun -AMAT151 --ntasks=32 -u -c7 --gpus-per-node=8  --ntasks-per-gpu=1 --gpu-bind=closest ./rmg-gpu input
