cmake_minimum_required (VERSION 2.6) 
project (RMG C Fortran) 

# Needed for a custom cmake module for FFTW
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Executable defined by file that contains main()
add_executable(rmg RMG/Common/main.c)


# static libraries are created for the object files in each subdir. Duplicates
# handle circular references
target_link_libraries (rmg Finite_diff RMG_Common)
target_link_libraries (rmg RMG_Common Finite_diff)
target_link_libraries (rmg RMG_Common Force)
target_link_libraries (rmg Force RMG_Common)
target_link_libraries (rmg Input RMG_Input)
target_link_libraries (rmg MG)
target_link_libraries (rmg Misc)
target_link_libraries (rmg Finite_diff Misc)
target_link_libraries (rmg Misc Finite_diff)
target_link_libraries (rmg RMG_Input Input)
target_link_libraries (rmg XC)
target_link_libraries (rmg RMG_Spin XC)
target_link_libraries (rmg XC RMG_Spin)
target_link_libraries (rmg RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg US_PP)
target_link_libraries (rmg RMG_Common)
target_link_libraries (rmg Force)
target_link_libraries (rmg Force US_PP)
target_link_libraries (rmg US_PP Force)


# Search for various libraries that we need
set(MPI_C_COMPILER "/opt/openmpi/bin/mpicc")
find_package(MPI REQUIRED)
if(MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_DIRS})
  target_link_libraries (rmg ${MPI_C_LIBRARIES})
endif (MPI_C_FOUND)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  target_link_libraries (rmg ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

# These are hard coded for now
target_link_libraries (rmg "/usr/local/lib64/scalapack/libscalapack.a")
target_link_libraries (rmg "/mnt/ssd2/briggs/Projects/RMG_RELEASE1/codes/lib/libxc/lib64/libxc.a")
target_link_libraries (rmg "/usr/local/atlas/lib/libtatlas.so")
target_link_libraries (rmg "/usr/lib64/libgfortran.so.3")

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  target_link_libraries (rmg ${BLAS_LIBRARIES})
endif (BLAS_FOUND)

if (${UNIX})
    set (CMAKE_THREAD_PREFER_PTHREAD)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
      target_link_libraries (rmg ${CMAKE_THREAD_LIBS_INIT})
    endif(CMAKE_USE_PTHREADS_INIT)
endif(${UNIX})

#find_package(FFTW REQUIRED)
#include_directories(${FFTW_INCLUDES})
#target_link_libraries (rmg ${FFTW_LIBRARIES})
target_link_libraries (rmg "/usr/local/lib/libfftw.a")

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall -O -g -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -mavx -fopenmp")        ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c")


# include files
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
include_directories("${RMG_BINARY_DIR}/RMG/Headers/")
include_directories("/opt/openmpi/include")


# Code subdirectories
add_subdirectory (Finite_diff) 
add_subdirectory (Force) 
add_subdirectory (Input) 
add_subdirectory (MG) 
add_subdirectory (Misc) 
add_subdirectory (RMG/Common) 
add_subdirectory (RMG/Input) 
add_subdirectory (RMG/Spin) 
add_subdirectory (RMG/Spin/XC) 
add_subdirectory (US_PP) 
add_subdirectory (XC) 

# Link subdirectories
link_directories (${RMG_BINARY_DIR}/Finite_diff) 
link_directories (${RMG_BINARY_DIR}/Force)
link_directories (${RMG_BINARY_DIR}/Input)
link_directories (${RMG_BINARY_DIR}/MG)
link_directories (${RMG_BINARY_DIR}/Misc)
link_directories (${RMG_BINARY_DIR}/RMG/Common)
link_directories (${RMG_BINARY_DIR}/RMG/Input)
link_directories (${RMG_BINARY_DIR}/RMG/Spin)
link_directories (${RMG_BINARY_DIR}/RMG/Spin/XC)
link_directories (${RMG_BINARY_DIR}/US_PP)
link_directories (${RMG_BINARY_DIR}/XC)


# configure make_conf.h to pass some of the CMake settings
# to the source code
option (RMG_EXPERIMENTAL_FEATURES "Turn ON to enable experimental features." OFF)
if (RMG_EXPERIMENTAL_FEATURES)
    unset (RMG_FAST_ORTHO CACHE)
    set (EXPERIMENTAL_FEATURES 1)
else (RMG_EXPERIMENTAL_FEATURES)
    unset (RMG_FAST_ORTHO CACHE)
    set (EXPERIMENTAL_FEATURES 0)
endif(RMG_EXPERIMENTAL_FEATURES)

option (RMG_HYBRID_MODEL "To enable MPI/PThreads hybrid model." ON)
if (RMG_HYBRID_MODEL)
    set (HYBRID_MODEL 1)
else (RMG_HYBRID_MODEL)
    set (HYBRID_MODEL 0)
endif(RMG_HYBRID_MODEL)

option (RMG_GAMMA_PT "Turn on for a gamma point only calculation." ON)
if (RMG_GAMMA_PT)
    set (GAMMA_PT 1)
else (RMG_GAMMA_PT)
    set (GAMMA_PT 0)
endif(RMG_GAMMA_PT)

option (RMG_FDIFF_BETA 
"Set this to 1 if you want to use finite difference method for calculating
 derivatives of beta. This is faster since it avoids doing 3 backwards fourier
 transforms per ion, but it may not be very accurate since the finite diff
 derivative is done on the coarse grid.
 Leave this set to 0 unless you know what you are doing " OFF)
if (RMG_FDIFF_BETA)
    set (FDIFF_BETA 1)
else (RMG_FDIFF_BETA)
    set (FDIFF_BETA 0)
endif(RMG_FDIFF_BETA)

option (RMG_MD_TIMERS "Extra fine timers, may cause a small slowdown but they are generally useful." ON)
if (RMG_MD_TIMERS)
    set (MD_TIMERS 1)
else (RMG_MD_TIMERS)
    set (MD_TIMERS 0)
endif(RMG_MD_TIMERS)

option (RMG_FAST_MEHR "Fast Mehrstellen operator." ON)
if (RMG_FAST_MEHR)
    set (FAST_MEHR 1)
else (RMG_FAST_MEHR)
    set (FAST_MEHR 0)
endif(RMG_FAST_MEHR)

set (FAST_ORTHO 0)
if(EXPERIMENTAL_FEATURES)
    option (RMG_FAST_ORTHO "Experimental fast ortho. Disabled by default." OFF)
    if (RMG_FAST_ORTHO)
        set (FAST_ORTHO 1)
    else (RMG_FAST_ORTHO)
        set (FAST_ORTHO 0)
    endif(RMG_FAST_ORTHO)
endif(EXPERIMENTAL_FEATURES)

option (RMG_GPU_ENABLED "Use GPU acceleration." OFF)
if (RMG_GPU_ENABLED)
    set (GPU_ENABLED 1)
else (RMG_GPU_ENABLED)
    set (GPU_ENABLED 0)
endif(RMG_GPU_ENABLED)

option (RMG_MAGMA_LIBS "Use MAGMA libraries." OFF)
if (RMG_MAGMA_LIBS)
    set (MAGMA_LIBS 1)
else (RMG_MAGMA_LIBS)
    set (MAGMA_LIBS 0)
endif(RMG_MAGMA_LIBS)

option (RMG_ASYNC_TRADES "Use asynchronous trade images." ON)
if (RMG_ASYNC_TRADES)
    set (ASYNC_TRADES 1)
else (RMG_ASYNC_TRADES)
    set (ASYNC_TRADES 0)
endif(RMG_ASYNC_TRADES)

# Send the configure file to the build directory
configure_file (
  "${PROJECT_SOURCE_DIR}/make_conf.h.in"
  "${PROJECT_BINARY_DIR}/RMG/Headers/make_conf.h"
  )


