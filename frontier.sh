#!/bin/bash
module load PrgEnv-gnu/8.6.0
module load gcc-native/13.2
module load cmake
module load cray-mpich/9.0.1
module load Core/24.00
module load bzip2
module load boost/1.85.0
module load craype-x86-milan
module load cray-fftw
module load cray-hdf5-parallel
module load rocm/6.2.4
export MPICH_GPU_SUPPORT_ENABLED=0

rm -rf build-frontier-gpu
mkdir build-frontier-gpu
cd build-frontier-gpu
export HIP_PATH=/opt/rocm-6.2.4/
export CXX=/opt/cray/pe/craype/2.7.33/bin/CC
cmake .. -DRMG_HIP_ENABLED=1 -DUSE_NCCL_OR_RCCL=1 -DUSE_ELPA_LIBS=1 -DCMAKE_PREFIX_PATH=/lustre/orion/mat739/proj-shared/elbriggs/usr/
make rmg-gpu -j 20 -k > rmg-gpu.log_1 2>&1
make rmg-gpu -j 20 -k > rmg-gpu.log_2 2>&1
make rmg-on-gpu -j 20 -k > rmg-on-gpu.log 2>&1
make rmg-negf-gpu -j 20 -k > rmg-negf-gpu.log 2>&1


