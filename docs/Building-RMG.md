RMG uses [Cmake](http://www.cmake.org) to manage the build process.[Cmake](http://www.cmake.org) is open source software and is available for a wide range of platforms. Building RMG is nearly automatic on most recent Linux distributions as long as the required packages are installed. Building on various cluster and supercomputer platforms often requires some additional work to get the environment properly configured. We will attempt to provide up to date instructions for various platforms at the the [issue tracker](https://github.com/RMGDFT/rmgdft/issues?q=is%3Aopen+is%3Aissue+label%3A%22Build+instructions%22).

### Supported platforms
* Most modern linux distributions on X86_64 hardware platforms.
* Power PC platforms.
* Cray Systems.
* Microsoft windows using the Windows Subsystem for Linux layer.

### Minimum prerequisites
* Cmake v3.1 (earlier versions may work but have not been extensively tested)
* Hardware platform with a 64 bit architecture.
* C++ compiler that supports the C++17 standard. GCC versions 7 and higher are known to work.
* Boost libraries. Versions 1.61 and higher.
* Blas libraries.
* Message passing (MPI) library that supports MPI_THREAD_SERIALIZED.

### Optional packages
* Cuda GPU computing libraries version 9.0 or higher.
* ROCM/HIP computing libraries version 5.4 or higher.
* Openbabel chemistry toolbox.
* PLplot scientific plotting package.

### Build instructions
Obtain the code from the RMG git [repository](https://github.com/RMGDFT/rmgdft)
If you have git installed on your local machine use <BR><BR>
`git clone https://github.com/RMGDFT/rmgdft.git`<BR><BR>
this will download the most recent code (not necessarily the most stable code).<BR><BR>
Alternatively you can download a tar or zip archive of a tagged release from [releases](https://github.com/RMGDFT/rmgdft/releases)
After downloading/unpacking the code or cloning it cd into the top level directory and execute the following instructions from the command line.<BR><BR>
`mkdir build`<BR>
`cd build`<BR>
`cmake ..`<BR><BR>
If you are trying to build a CUDA enabled version of RMG change the last line to<BR><BR>
`cmake -DRMG_CUDA_ENABLED=1 ..`<BR><BR>
If you are trying to build a ROCM/HIP enabled version of RMG change the last line to<BR><BR>
`cmake -DRMG_HIP_ENABLED=1 ..`<BR><BR>
If cmake completes without errors you can now try building RMG using <BR>                             
`make -jN rmg-cpu`<BR><BR>
or if you have used -DRMG_CUDA_ENABLED=1 or -DRMG_HIP_ENABLED=1<BR><BR>
`make -jN rmg-gpu`<BR><BR>
where N is the number of CPU cores to use in a parallel build.

### Debug builds
To generate a debug build edit the main CMakeList.txt file in the root directory of the distribution and change<BR><BR>
` set(CMAKE_BUILD_TYPE Release)`<BR>
`#set(CMAKE_BUILD_TYPE Debug)`<BR>
to<BR>
`#set(CMAKE_BUILD_TYPE Release)`<BR>
`set(CMAKE_BUILD_TYPE Debug)`<BR>

### Testing
After you have built the rmg-cpu executable you can run system level tests by executing.<BR>
`make test`<BR>
In order to successfully complete all of the tests your machine should have at least 16 GBytes of RAM. The test script will try to use as many CPU cores as your test machine has up to a maximum of 24. On a 24 core Ryzen threadripper it takes around 8 minutes to complete all of the tests.
 
### Troubleshooting
If your platform supports modules you may not have the correct set of modules loaded. Try the command `module avail` to see if your system supports modules and what is available.
If you run into problems you can ask for help at the [issue tracker.](https://github.com/RMGDFT/rmgdft/issues?q=is%3Aopen+is%3Aissue+label%3A%22Build+instructions%22)<BR>

### System specific: Frontier at ORNL
As of 07/09/23 the following module and cmake command works on Frontier<BR><BR>
`module load cmake`<BR>
`module load PrgEnv-gnu/8.3.3`<BR>
`module load bzip2`<BR>
`module load boost/1.79.0`<BR>
`module load cray-fftw`<BR>
`module load cray-hdf5-parallel`<BR>
`module load craype-accel-amd-gfx90a`<BR>
`module load rocm/5.4.3`<BR>
`export MPICH_GPU_SUPPORT_ENABLED=0`<BR>
`cmake -DRMG_HIP_ENABLED=1 ..`<BR>
`make -j32 rmg-gpu`<BR><BR>
As of 02/13/24 the following may be used to build the develop branch of RMG using ROCM 6.0 on Frontier<BR><BR>
`module load cmake`<BR>
`module load PrgEnv-gnu-amd/8.5.0`<BR>
`module load bzip2`<BR>
`module load boost/1.79.0`<BR>
`module load cray-fftw`<BR>
`module load cray-hdf5-parallel`<BR>
`module load craype-accel-amd-gfx90a`<BR>
`module load amd-mixed/6.0.0`<BR>
`export MPICH_GPU_SUPPORT_ENABLED=0`<BR>
`cmake -DRMG_HIP_ENABLED=1 ..`<BR>
`make -j32 rmg-gpu`<BR>

### System specific: Summit at ORNL
As of 7/11/23 the following environment, module and cmake command work on Summit<BR> 
`export FC=/sw/summit/gcc/9.1.0-alpha+20190716/bin/gfortran`<BR>
`export CC=/sw/summit/gcc/9.1.0-alpha+20190716/bin/gcc`<BR>
`export CXX=/sw/summit/gcc/9.1.0-alpha+20190716/bin/g++`<BR>
`export BLA_VENDOR=OpenBLAS`<BR>
`module load gcc`<BR>
`module load boost`<BR>
`module load openblas/0.3.15-omp`<BR>
`module load cuda`<BR>
`module load fftw`<BR>
`module load hdf5`<BR>
`module load cmake/3.20.2`<BR>
`module load bzip2`<BR>
`cmake -DRMG_CUDA_ENABLED=1 -DUSE_INTERNAL_SCALAPACK=1 ..`<BR>
`make -j32 rmg-gpu`<BR>

### System specific: NERSC machine Perlmutter
As of 07/09/23 the following module and cmake command works on Perlmutter<BR><BR>
`module load cmake`<BR>
`module load cray-fftw`<BR>
`module load cray-hdf5-parallel`<BR>
`cmake -DRMG_CUDA_ENABLED=1 -DCMAKE_Cuda_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/bin/nvcc ..`<BR>
`make -j32 rmg-gpu`