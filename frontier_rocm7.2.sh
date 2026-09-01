#!/bin/bash
# frontier build script for rocm 7.2.0 and ELPA
module reset
module load PrgEnv-gnu/8.7.0   
module load cpe/26.03
module load cray-mpich/9.1.0
module load rocm/7.0.2
module swap rocm/7.0.2 rocm/7.2.0
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load boost/1.88.0
module load cray-fftw
module load cray-hdf5-parallel
export CC=cc
export CXX=CC
export FC=ftn
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export MPICH_GPU_SUPPORT_ENABLED=0

rm -rf build-frontier-gpu
mkdir build-frontier-gpu
cd build-frontier-gpu
cmake -DRMG_HIP_ENABLED=1 -DHIP_PATH=/opt/rocm-7.2.0 -DCRAY_HOST=1 -DUSE_ELPA_LIBS=1 -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH;/lustre/orion/mat739/proj-shared/elbriggs/elpa_rocm7.2/" ..

make rmg-gpu -j 20 -k > rmg-gpu.log_1 2>&1
make rmg-gpu -j 20 -k > rmg-gpu.log_2 2>&1
make rmg-on-gpu -j 20 -k > rmg-on-gpu.log 2>&1
make rmg-negf-gpu -j 20 -k > rmg-negf-gpu.log 2>&1


