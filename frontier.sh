#!/bin/bash
module load PrgEnv-gnu/8.6.0
module load gcc-native/13.2
module load cmake
module load Core/24.00
module load bzip2
module load boost/1.85.0
module load craype-x86-milan
module load cray-fftw
module load cray-hdf5-parallel
module load craype-accel-amd-gfx90a
module load rocm/6.3.1
export MPICH_GPU_SUPPORT_ENABLED=0

rm -rf build-frontier-gpu
mkdir build-frontier-gpu
cd build-frontier-gpu
cmake .. -DRMG_HIP_ENABLED=1 -DHIP_PATH="/opt/rocm-6.3.1/"
make rmg-gpu -j 20 -k > rmg-gpu.log 2>&1
make rmg-on-gpu -j 20 -k > rmg-on-gpu.log 2>&1
make rmg-negf-gpu -j 20 -k > rmg-negf-gpu.log 2>&1


