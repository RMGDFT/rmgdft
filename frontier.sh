#!/bin/bash
module load PrgEnv-gnu/8.5.0
module load cmake
module load boost/1.79.0
module load craype-x86-milan
module load cray-fftw
module load cray-hdf5-parallel
module load craype-accel-amd-gfx90a
module load rocm/6.0.0
module load libfabric/1.15.2.0
export MPICH_GPU_SUPPORT_ENABLED=0
export FI_MR_CACHE_MONITOR=disabled

rm -rf build-frontier-gpu
mkdir build-frontier-gpu
cd build-frontier-gpu
cmake .. -DRMG_HIP_ENABLED=1 -DHIP_PATH="/opt/rocm-6.0.0/"
make rmg-gpu -j 20
make rmg-on-gpu -j 20
make rmg-negf-gpu -j 20


