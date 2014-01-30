#
#    RMG master cmake file
#
#    Subprojects:
#        rmg        Base code
#        rmg-on     ON code
#        NEGF
#
# Environment variables that may also control configuration
# MKLROOT = root directory of Intel MKL libraries if installed.
#
#
#
cmake_minimum_required (VERSION 2.6) 
project (RMG C Fortran) 


# Executables defined by the file that contains main()
add_executable(rmg RMG/Common/main.c)		# Base code
add_executable(rmg-on ON/Common/md.c)		# Order N


# Require a 64 bit system
if (NOT CMAKE_SIZEOF_VOID_P MATCHES 8 )
    message( SEND_ERROR "RMG requires a 64 bit architechture." )
endif (NOT CMAKE_SIZEOF_VOID_P MATCHES 8)


# Require pthreads on UNIX
if (${UNIX})
    set (CMAKE_THREAD_PREFER_PTHREAD)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
        target_link_libraries (rmg ${CMAKE_THREAD_LIBS_INIT})
    endif(CMAKE_USE_PTHREADS_INIT)
endif(${UNIX})

# Custom cmake module for FFTW
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)




# static libraries are created for the object files in each subdir. Duplicates
# are needed to handle circular references
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



# Search for an MPI implementation
#set(MPI_C_COMPILER "/opt/openmpi/bin/mpicc")
find_package(MPI REQUIRED)
if(MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  target_link_libraries (rmg ${MPI_C_LIBRARIES})
endif (MPI_C_FOUND)


# For various math libraries start by looking for Intel MKL. No good cmake module currently for MKL.
# Check if environment variable is set
set(MKLROOT $ENV{MKLROOT})
if (NOT MKLROOT)
    # Check the default location
    if (EXISTS "/opt/intel/mkl")
        set(MKLROOT "/opt/intel/mkl")
    endif (EXISTS "/opt/intel/mkl")
endif (NOT MKLROOT)

if(MKLROOT)
    set(MKL_LIBDIR "${MKLROOT}/lib/intel64")
    set(MKL_LIB64 "${MKLROOT}/lib/intel64")
    if(EXISTS "${MKL_LIB64}")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_scalapack_lp64.a" "${MKL_LIB64}/libmkl_intel_lp64.a")
#        target_link_libraries (rmg "/usr/local/intel_fftw/libfftw2xc_double_gnu.a")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_cdft_core.a")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_intel_lp64.a")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_core.a")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_gnu_thread.a" "${MKL_LIB64}/libmkl_intel_lp64.a")
        target_link_libraries (rmg "${MKL_LIB64}/libmkl_blacs_openmpi_lp64.a")
        set(RMG_MATH_LIBS "done")
    else(EXISTS "${MKL_LIB64}")
        message( SEND_ERROR "RMG requires 64 bit MKL." )
    endif(EXISTS "${MKL_LIB64}")
endif(MKLROOT)


# MKL includes blas, blacs, lapack and scalapack. If it is not present we have to look for these libraries
# individually
if(NOT RMG_MATH_LIBS)
    find_package(BLAS REQUIRED)
    if(BLAS_FOUND)
        target_link_libraries (rmg ${BLAS_LIBRARIES})
    endif (BLAS_FOUND)

    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND)
        target_link_libraries (rmg ${LAPACK_LIBRARIES})
    endif (LAPACK_FOUND)
endif(NOT RMG_MATH_LIBS)


# Look for GPU libraries. Only CUDA supported for now
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
  target_link_libraries(rmg ${CUDA_CUDART_LIBRARY})
  target_link_libraries(rmg ${CUDA_CUDA_LIBRARY})
  target_link_libraries(rmg ${CUDA_cublas_LIBRARY})
  include_directories(${CUDA_INCLUDE_DIRS})
endif(CUDA_FOUND)

# These are hard coded for now
target_link_libraries (rmg "/home/emil/Projects/RMG_RELEASE1/codes/lib/libxc/lib64/libxc.a")

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
target_link_libraries (rmg ${FFTW_LIBRARIES})


# Compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall -O3 -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -mavx -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
    set(CMAKE_EXE_LINKER_FLAGS "-ldl")  ## Strip binary
endif()
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c")


# include files
include_directories("${RMG_SOURCE_DIR}/Headers/")
#include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
#include_directories("${RMG_BINARY_DIR}/RMG/Headers/")
#include_directories("/opt/openmpi/include")


# Common code subdirectories
add_subdirectory (Finite_diff) 
add_subdirectory (Force) 
add_subdirectory (Input) 
add_subdirectory (MG) 
add_subdirectory (Misc) 
add_subdirectory (US_PP) 
add_subdirectory (XC) 

# Target specific code subdirectories
add_subdirectory (RMG/Common EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Input EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin/XC EXCLUDE_FROM_ALL) 

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
    if(CUDA_FOUND)
        set (GPU_ENABLED 1)
    else(CUDA_FOUND)
        set (GPU_ENABLED 0)
    endif(CUDA_FOUND)
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
  "${PROJECT_BINARY_DIR}/Headers/make_conf.h"
  )


