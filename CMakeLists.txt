#
#    RMG master cmake file
#
#    Subprojects
#        rmg-cpu        Base code cpu only
#        rmg-gpu        Base code with gpu support
#        rmg-on-cpu     ON code cpu only
#        rmg-on-gpu     ON code with gpu support
#        rmg-negf-cpu   Non equilibrium greens function code cpu only
#        rmg-negf-gpu   Non equilibrium greens function code with gpu support

#
# Environment variables that may also control configuration
#
# MKLROOT = root directory of Intel MKL libraries if installed.
# BOOST_ROOT = root directory of Boost installation
# CC  = Full path to C compiler
# CXX = Full path to C++ compiler
# FC  = Full path to Fortran compiler
#
# On some platforms an environment variable may need to be set to
# find an optimized BLAS library. For example on Summit at ORNL
# with the cmake/3.14.2 module loaded one needs to set
# export BLA_VENDOR=IBMESSL in order for cmake to find the essl libs.
#


cmake_minimum_required (VERSION 3.1) 


# force an out-of-source build
string( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" RMG_COMPILE_INPLACE )
if( RMG_COMPILE_INPLACE )
        message( FATAL_ERROR "Building RMG with CMake requires an out-of-source build. To proceed:
    rm -rf CMakeCache.txt CMakeFiles/   # delete files in ${CMAKE_SOURCE_DIR}
    mkdir build
    cd build
    cmake ..
    make -jX target
    Where -jX with X=ncpus to compile with and target may be rmg-cpu or rmg-gpu")
endif()


enable_language(Fortran)
project (RMG CXX C Fortran) 


# Extract revision information to include in executables
find_package(Git)
if(Git_FOUND)
    set(VERSION_COMMAND "git")
    set(VERSION_ARGS "describe")
    unset(git_version CACHE)
    execute_process (COMMAND "${VERSION_COMMAND}" "${VERSION_ARGS}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} RESULT_VARIABLE git_result OUTPUT_VARIABLE git_version)
    if(git_version)
        string(REPLACE "\n" ";" git_version1 ${git_version})
        string(REPLACE "-" "_" git_version_clean ${git_version1})
        add_definitions(-DRMG_REVISION="${git_version_clean}")
    else()
        add_definitions(-DRMG_REVISION="Unknown")
        set(git_version_clean "Unknown")
    endif()
else()
    add_definitions(-DRMG_REVISION="Unknown")
    set(git_version_clean "Unknown")
endif()


# Require a 64 bit system
if (NOT CMAKE_SIZEOF_VOID_P MATCHES 8 )
    message( SEND_ERROR "RMG requires a 64 bit architechture." )
endif()


# Check if this is a cray
if(DEFINED ENV{CRAY_LIBSCI_VERSION})
    # Cray platforms
    set (CRAY_HOST 1)
endif()

# User settable options from cmake-gui, ccmake for text terminals
# or command line
set(RMG_CUDA_ENABLED OFF CACHE BOOL "Enable CUDA support.")
set(RMG_HIP_ENABLED OFF CACHE BOOL "Enable HIP support.")
set(RMG_MAGMA_ENABLED OFF CACHE BOOL "Enable MAGMA support.")
set(QMCPACK_SUPPORT ON CACHE BOOL "Generates a QMCPACK restart file.")
set(USE_INTERNAL_PSEUDOPOTENTIALS ON CACHE BOOL "Include internal pseudopotentials in binary.")
set(USE_OPENBABEL_LIBS OFF CACHE BOOL "Openbabel support for importing atomic data.")
set(USE_PLPLOT_LIBS OFF CACHE BOOL "Use plplot for runtime graphs.")
set(USE_INTERNAL_SCALAPACK OFF CACHE BOOL "Use internal scalapack.")

# For profiling blas calls
option(BLAS_PROFILING "Enable blas profiling" OFF)
if(BLAS_PROFILING)
    add_definitions(-DBLAS_PROFILE=1)
else()
    add_definitions(-DBLAS_PROFILE=0)
endif()

# Build 12th order laplacian for testing purposes
option(TWELTH_ORDER_FD "Build 12th order laplacian for testing" OFF)
if(TWELFTH_ORDER_FD)
    add_definitions(-DTWELFTH_ORDER_FD=1)
endif()

if(RMG_CUDA_ENABLED AND RMG_HIP_ENABLED)
    message(FATAL_ERROR "\nRMG_CUDA_ENABLED and RMG_HIP_ENABLED are mutually exclusive. Pick one or the other for a GPU build. If neither is set a CPU only build will be performed.")
endif()


# Now that the scalapack source has been included in the RMG distribution this should
# always be available unless there is some issue with building it so we define it here.
add_definitions(-DSCALAPACK_LIBS=1)


# Helps to find the correct version of boost on some systems. In particular
# many supercomputers will have a system version located in the standard
# directories that does not match the version required by the compiler you
# want to use. If cmake is not picking up the correct version try toggling this.
option(Boost_NO_SYSTEM_PATHS "If you have multiple versions of boost installed toggling this may be required" OFF)
# If you are having problems with boost at link time try toggling these (e.g. some IBM Power machines)
option(Boost_USE_STATIC_LIBS        OFF)
option(Boost_USE_MULTITHREADED      ON)
option(USE_ELPA_LIBS "Experimental." OFF)
option(USE_ELEMENTAL_LIBS "Experimental." OFF)


# RMG requires c++14 and C99
include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)


# Check for c++14 and c99 support
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    set(CMAKE_CXX_STANDARD 14)
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set(CMAKE_CXX_STANDARD 14)
else()
    message( WARNING "The detected compiler ${CMAKE_CXX_COMPILER} doesn't support the c++14 standard or requires a flag cmake does not understand. RMG will probably not build correctly.")
endif()


CHECK_C_COMPILER_FLAG("-std=c99" COMPILER_SUPPORTS_C99)
CHECK_C_COMPILER_FLAG("-c99" COMPILER_SUPPORTS_DASH_C99)
if(COMPILER_SUPPORTS_C99)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_C_STANDARD 99)
elseif(COMPILER_SUPPORTS_DASH_C99)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")
    set(CMAKE_C_STANDARD 99)
else()
    message( WARNING "The detected compiler ${CMAKE_C_COMPILER} does not support the c99 standard or requires a flag cmake does not understand. RMG will probably not build correctly.")
endif()


# Standard modules for platform checks
include (CheckFunctionExists)
include (CheckSymbolExists)
include (CheckCSourceCompiles)
include (CheckCXXSourceCompiles)
include (CheckIncludeFiles)
include (CheckLibraryExists)
include (CMakeForceCompiler)
include (FortranCInterface)


# Fortran C interface mangling
FortranCInterface_VERIFY(CXX QUIET)
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/Headers/rmg_mangling_in.h
                         MACRO_NAMESPACE RMG_FC_
                         SYMBOL_NAMESPACE RMG_FC_)
# Used in Headers/rmg_mangling.h
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Cray)$")
    add_definitions(-DRMG_CRAY_FORTRAN=1)
endif()


# Modify this if you need some special compiler flags.
set(RMG_EXTRA_COMPILER_FLAGS "-Wno-deprecated-declarations")


# Generic Release or Debug
# This should normally be set automatically to release. It can also be set when cmake is invoked
# by using -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()


# Load custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)


# Search for an MPI implementation. On Crays cmake does not pick up some
# vars but the compiler wrappers take care of it automatically.
if(CRAY_HOST)
    find_package(MPI)
else()
    find_package(MPI REQUIRED)
endif()
set(RMGLIBS "${MPI_CXX_LIBRARIES}")
list(APPEND RMGLIBS ${MPI_C_LIBRARIES})


if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
if(MPI_C_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
endif()


# Look for scalapack
if(NOT USE_INTERNAL_SCALAPACK)
   if(NOT CRAY_HOST)
       find_package(SCALAPACK)
   endif()
endif()


# Get Openmp flags
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")


if(USE_INTERNAL_PSEUDOPOTENTIALS)
    add_definitions(-DINTERNAL_PP=1)
    list(APPEND RMGLIBS InternalPseudo)
endif()


if(RMG_CUDA_ENABLED)
    add_definitions(-DCUDA_ENABLED=1)
    set(rmg-negf "rmg-negf-gpu")
    set(rmg "rmg-gpu")
    set(rmg-on "rmg-on-gpu")
else()
    if(RMG_HIP_ENABLED)
        add_definitions(-DHIP_ENABLED=1)
        set(rmg-negf "rmg-negf-gpu")
        set(rmg "rmg-gpu")
        set(rmg-on "rmg-on-gpu")
    else()
        add_definitions(-DCUDA_ENABLED=0)
        add_definitions(-DMAGMA_LIBS=0)
        set(rmg-negf "rmg-negf-cpu")
        set(rmg "rmg-cpu")
        set(rmg-on "rmg-on-cpu")
    endif()
endif()

if(RMG_HIP_ENABLED)
    if(NOT DEFINED HIP_PATH)
        if(NOT DEFINED ENV{HIP_PATH})
            set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
        else()
            set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
        endif()
    endif()
    list (APPEND CMAKE_MODULE_PATH "${HIP_PATH}/cmake")
    find_package(HIP)
    if(HIP_FOUND)
        message(STATUS "Found HIP: " ${HIP_VERSION})
        add_definitions(-DHIP_ENABLED=1)
        include_directories(${HIP_INCLUDE_DIR})
        include_directories("${HIP_PATH}/include")
        if(DEFINED ENV{ROCM_PATH})
            set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
        else()
            set(ROCM_PATH "/opt/rocm")
        endif()
        include_directories("${ROCM_PATH}/hsa/include")
        include_directories("${ROCM_PATH}")
        include_directories("${ROCM_PATH}/include")
        if(DEFINED HIP_LIBRARY)
            list(APPEND RMGLIBS ${HIP_LIBRARY})
        else()
            list(APPEND RMGLIBS "${HIP_PATH}/lib/libamdhip64.so")
        endif()
    endif()

    find_package(HIPBLAS HINTS $ENV{OLCF_HIPBLAS_ROOT} $ENV{HIPBLAS_PATH} REQUIRED)
    if(HIPBLAS_FOUND)
        message(STATUS "Found HIPBLAS: " ${HIPBLAS_VERSION})
        list(APPEND RMGLIBS ${HIPBLAS_LIBRARY})
        include_directories("${HIPBLAS_INCLUDE_DIR}")
    endif()

    find_package(rocfft HINTS $ENV{OLCF_ROCFFT_ROOT} $ENV{ROCFFT_PATH} REQUIRED)
    if(rocfft_FOUND)
        message(STATUS "Found rocfft: " ${rocfft_VERSION})
        list(APPEND RMGLIBS ${rocfft_LIBRARY})
        include_directories("${rocfft_INCLUDE_DIR}")
    endif()

    find_package(rocsolver REQUIRED)
    if(rocsolver_FOUND)
        message(STATUS "Found rocsolver: " ${rocsolver_VERSION})
        list(APPEND RMGLIBS ${rocsolver_LIBRARY})
        include_directories("${rocsolver_INCLUDE_DIR}")
    endif()

    find_package(hipsolver REQUIRED)
    if(hipsolver_FOUND)
        message(STATUS "Found hipsolver: " ${hipsolver_VERSION})
        list(APPEND RMGLIBS ${hipsolver_LIBRARY})
        include_directories("${hipsolver_INCLUDE_DIR}")
    endif()

    if(NOT DEFINED __HIP_PLATFORM_HCC__)
        add_definitions(-D__HIP_PLATFORM_HCC__=1)
    endif()
    if(RMG_MAGMA_ENABLED)
        find_package(MAGMA)
        if(MAGMA_FOUND)
            add_definitions(-DMAGMA_LIBS=1)
            include_directories(${MAGMA_INCLUDES})
            list(APPEND RMGLIBS ${MAGMA_LIBRARIES})
        endif()
    endif()
endif()

# Debug prints out all vars
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# System specific flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(OperatingSystem "Linux")
    add_definitions(-DLINUX=1)
    set(ZLIB "z")
    set(RTLIB "rt")
    list(APPEND RMGLIBS z)
    list(APPEND RMGLIBS rt)
    list(APPEND RMGLIBS m)
endif()


# GNUCC specific compiler flags
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pthread -m64 -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -Wno-write-strings -pthread -m64 -fext-numeric-literals -fopenmp")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
        set (CMAKE_Fortran_FLAGS "-lg2c -cpp -c -fopenmp -fallow-argument-mismatch")
        set (CMAKE_Fortran_FLAGS_RELEASE "-lg2c -cpp -c -fopenmp -fallow-argument-mismatch")
    else()
        set (CMAKE_Fortran_FLAGS "-lg2c -cpp -c -fopenmp")
        set (CMAKE_Fortran_FLAGS_RELEASE "-lg2c -cpp -c -fopenmp")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS "-ldl -fopenmp")
    set(RMG_EXTRA_COMPILER_FLAGS "${RMG_EXTRA_COMPILER_FLAGS} -fcx-limited-range -fno-trapping-math")
    list(APPEND RMGLIBS z)
    list(APPEND RMGLIBS rt)
    list(APPEND RMGLIBS m)
endif()


# MKL specific compiler flags
#if(CMAKE_C_COMPILER_ID)
#    # Intel compiler
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX -D_NOMINMAX -I${MKLROOT}/include")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNOMINMAX -D_NOMINMAX -I${MKLROOT}/include")
#    #set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")
#endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${RMG_EXTRA_COMPILER_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RMG_EXTRA_COMPILER_FLAGS}")


# This file contains system specific paths. If cmake is unable to find
# a specific package or library look in this file and uncomment or add
# specific search paths.
include(CMakeFindRootPath.inc)


# adds qmcpack support to base code
if(QMCPACK_SUPPORT)
    add_definitions(-DQMCPACK_SUPPORT=1)
endif()


# Require pthreads on UNIX
if (${UNIX})
    find_package(Threads REQUIRED)
    list(APPEND RMGLIBS ${RTLIB})
    list(APPEND RMGLIBS ${CMAKE_THREAD_LIBS_INIT})
endif()


find_package(NUMA)
if(NUMA_FOUND AND NUMA_INCLUDES)
    include_directories(${NUMA_INCLUDES})
    add_definitions(-DUSE_NUMA=1)
    list(APPEND RMGLIBS ${NUMA_LIBRARIES})
endif()


find_package(HWLOC)
if(HWLOC_FOUND AND HWLOC_INCLUDES)
    include_directories(${HWLOC_INCLUDES})
    #add_definitions(-DUSE_HWLOC=1)
    list(APPEND RMGLIBS ${HWLOC_LIBRARIES})
endif()


#experimental
if(USE_ELPA_LIBS)
    message(WARNING "ELPA is experimental and requires manual editing of CMakeLists.txt to enable.")
    #add_definitions(-DUSE_ELPA=1)
    #list(APPEND RMGLIBS /home/briggs/src/elpa-2016.11.001.pre/.libs/libelpa.so)
    #list(APPEND RMGLIBS /usr/lib64/libelpa.so)
endif()


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "^(FreeBSD|DragonFly)$")
    if(CRAY_HOST)
        add_definitions(-DHAVE_ASYNC_ALLREDUCE=1)
    else(CRAY_HOST)
        include("PlatformChecks/CheckAsyncAllreduce.inc")
        if(RMG_ASYNC_ALLREDUCE_INTERNAL)
            add_definitions(-DHAVE_ASYNC_ALLREDUCE=1)
        else()
            add_definitions(-DHAVE_ASYNC_ALLREDUCE=0)
        endif()
    endif()
endif()


if(USE_PLPLOT_LIBS)
    find_package(PLPLOT)
    if(PLPLOT_FOUND)
        include("PlatformChecks/CheckPlplot.inc")
        if(USE_PLPLOT_LIBS_INTERNAL)
            include_directories(${PLplot_INCLUDE_DIR})
            include_directories("${PLplot_INCLUDE_DIR}/plplot")
            add_definitions(-DPLPLOT_LIBS=1)
            list(APPEND RMGLIBS ${PLplot_cxx_LIBRARY})
            list(APPEND RMGLIBS ${PLplot_LIBRARY})
        else()
            add_definitions(-DPLPLOT_LIBS=0)
        endif()
    else()
        add_definitions(-DPLPLOT_LIBS=0)
    endif()
endif()


set(Boost_NO_BOOST_CMAKE on)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(BOOST_ROOT "...")
    if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT $ENV{BOOST_ROOT})
    endif()
    find_package(Boost 1.61 REQUIRED COMPONENTS thread system iostreams program_options filesystem REQUIRED)
endif()
list(APPEND RMGLIBS ${Boost_LIBRARIES})
add_definitions(-DBOOST_FILESYSTEM_VERSION=3)


if(CRAY_HOST)
    find_package(BZ2 REQUIRED)
    find_package(Z REQUIRED)
    list(APPEND RMGLIBS  ${BZ2_LIBRARIES})
    list(APPEND RMGLIBS  ${Z_LIBRARIES})
endif()


# Look for openbabel
if(USE_OPENBABEL_LIBS)
    find_package(OPENBABEL)
    include("PlatformChecks/CheckOpenbabel.inc")

    if(CHECK_OPENBABEL_INTERNAL)
        include_directories(${OPENBABEL_INCLUDES})
        add_definitions(-DOPENBABEL_LIBS=1)
        list(APPEND RMGLIBS ${OPENBABEL_LIBRARIES})
    else()
        add_definitions(-DOPENBABEL_LIBS=0)
    endif()
endif()


# For various math libraries start by looking for Intel MKL. No good cmake module currently for MKL.
# Check if environment variable is set. Unfortunately on many clusters there are multiple versions
# of MKL installed and this does not always work correctly.
#if(NOT MKLROOT)
#    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#        # Check the default location
#        if (EXISTS "/opt/intel/mkl")
#            set(MKLROOT "/opt/intel/mkl")
#        endif()
#    endif()
#endif()


if(NOT CRAY_HOST)
    if(DEFINED ENV{MKLROOT})
	set(MKLROOT $ENV{MKLROOT})
        include_directories(${MKLROOT}/include)
        add_definitions(-DUSE_MKL)
        set(FFTW_INCLUDES "${MKLROOT}/include/fftw")
        set(MKL_LIBDIR "${MKLROOT}/lib/intel64")
        set(MKL_LIB64 "${MKLROOT}/lib/intel64")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            if(EXISTS "${MKL_LIB64}")
                list(APPEND RMGLIBS "${MKL_LIB64}/libmkl_intel_lp64.a")
                list(APPEND RMGLIBS "${MKL_LIB64}/libmkl_core.a")
                list(APPEND RMGLIBS "${MKL_LIB64}/libmkl_gnu_thread.a")
                list(APPEND RMGLIBS "${MKLROOT}/lib/intel64/libmkl_intel_lp64.a")
                list(APPEND RMGLIBS "${MKLROOT}/lib/intel64/libmkl_gnu_thread.a")
                list(APPEND RMGLIBS "${MKLROOT}/lib/intel64/libmkl_core.a")
                set(RMG_MATH_LIBS "done")
                set(RMG_MKLLIBS 1)
            else()
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif()
        endif()
    endif()
    if(USE_INTERNAL_SCALAPACK)
        list(APPEND RMGLIBS scalapack)
    else()
	if(SCALAPACK_FOUND)
            list(APPEND RMGLIBS ${SCALAPACK_LIBRARIES})
	else()
            list(APPEND RMGLIBS scalapack)
        endif()
    endif()
endif()


# MKL includes fftw blas, blacs, lapack and scalapack. If it is not present we have to look for these libraries
# individually
if(NOT ${FFTW_INCLUDES})
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDES})
    list(APPEND RMGLIBS ${FFTW_LIBRARIES})
    list(APPEND RMGLIBS ${FFTWF_LIBRARIES})
endif()


if(NOT RMG_MATH_LIBS)
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDES})
    list(APPEND RMGLIBS ${FFTW_LIBRARIES})
    list(APPEND RMGLIBS ${FFTWF_LIBRARIES})

    if(NOT $ENV{BLA_VENDOR} STREQUAL "")
        set(BLA_VENDOR $ENV{BLA_VENDOR})
    else()
        if(NOT BLA_VENDOR)
            if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64.*")
                message(STATUS "Trying to find Openblas")
                SET(BLA_VENDOR OpenBLAS)
            elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
                message(STATUS "Trying to find BLAS from ARM Performance Library")
                SET(BLA_VENDOR Arm_mp)
            endif() 
        endif()
    endif()

    if(NOT CRAY_HOST)
        find_package(BLAS REQUIRED)
        if(BLAS_FOUND)
            list(APPEND RMGLIBS ${BLAS_LIBRARIES})
        endif()
    endif()
endif()

# Check blas libraries
include("PlatformChecks/BlasFilter.inc")

# HDF5 libraries for QMCPACK integration
set(ENABLE_PHDF5 1 CACHE BOOL "Enable phdf5 for output")
set(HAVE_LIBHDF5 0)
set(HDF5_PREFER_PARALLEL 0)
if(HAVE_MPI)
   set(HDF5_PREFER_PARALLEL 1)
endif()
find_package(HDF5 COMPONENTS C)
if(HDF5_FOUND)
    set(HAVE_LIBHDF5 1)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    if ( CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG )
        if ( CMAKE_BUILD_TYPE MATCHES DEBUG )
            set( HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG} )
        else()
            set( HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE} )
        endif()
    endif()
    if(HDF5_IS_PARALLEL)
        if(ENABLE_PHDF5)
            message(STATUS "Using HDF5 parallel collective I/O")
        else()
            message(STATUS "Using HDF5 non-scalable serial I/O")
        endif()
    else()
        if(ENABLE_PHDF5)
            message(STATUS "Using HDF5 non-scalable serial I/O due to the lack of library support for parallel")
            set(ENABLE_PHDF5 0)
        else()
            message(STATUS "Using HDF5 non-scalable serial I/O")
        endif()
    endif()
    list(APPEND RMGLIBS ${HDF5_LIBRARIES})
else()
    message(FATAL_ERROR "Require hdf5 1.6.4 or higher. Set HDF5_ROOT")
endif()


# Look for GPU libraries. Only CUDA supported for now
if(RMG_CUDA_ENABLED)
    find_package(CUDA REQUIRED)
    message( WARNING "CUDA_VERSION_MAJOR = ${CUDA_VERSION_MAJOR}")
    message( WARNING "CUDA_VERSION_MINOR = ${CUDA_VERSION_MINOR}")
    if(${CUDA_VERSION_MAJOR} LESS 6)
        message(FATAL_ERROR "\nCompiling RMG with GPU support requires CUDA version 6.0 or greater. Either set CUDA_ENABLED=0 to build a CPU only version or upgrade your cuda installation.") 
    endif()

    if(NOT CRAY_HOST)
        list(APPEND RMGLIBS ${CUDA_CUDA_LIBRARY})
        list(APPEND RMGLIBS ${CUDA_CUDART_LIBRARY})
        list(APPEND RMGLIBS ${CUDA_cusolver_LIBRARY})
        list(APPEND RMGLIBS ${CUDA_cufft_LIBRARY})
    endif()
    list(APPEND RMGLIBS ${CUDA_cublas_LIBRARY})
    list(APPEND RMGLIBS ${CUDA_cusparse_LIBRARY})
    list(APPEND RMGLIBS ${CUDA_cusolver_LIBRARY})
    list(APPEND RMGLIBS ${CUDA_cufft_LIBRARY})
    include_directories(${CUDA_INCLUDE_DIRS})
    find_library (cusolverMg_LIBRARY NAMES cusolverMg PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
    list(APPEND RMGLIBS ${cusolverMg_LIBRARY})
    # Cuda functions that we build are in Gpufuncs/CMakeLists.txt
    add_definitions(-DCUDA_VERSION_MAJOR=${CUDA_VERSION_MAJOR})
    add_definitions(-DCUDA_VERSION_MINOR=${CUDA_VERSION_MINOR})
endif()


if(USE_ELEMENTAL_LIBS)
    add_definitions(-DELEMENTAL_LIBS=1)
    find_package(ELEMENTAL)
    if(NOT ELEMENTAL_FOUND)
        execute_process(
        COMMAND tar xvf ${PROJECT_SOURCE_DIR}/lib/Elemental-0.85.tgz -C ${PROJECT_SOURCE_DIR}/lib )
        execute_process(
        COMMAND mkdir ${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build)
        execute_process(
        COMMAND cmake .. 
            WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build)
        execute_process(COMMAND make  -j 16
            WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build)
    endif()

    find_package(ELEMENTAL REQUIRED)
    include_directories(${ELEMENTAL_INCLUDES})
    list(APPEND RMGLIBS ${ELEMENTAL_LIBRARIES})
    find_package(PMRRR)
    include_directories(${PMRRR_INCLUDES})
    list(APPEND RMGLIBS ${PMRRR_LIBRARIES})
    list(APPEND RMGLIBS ${PMRRR_LIBRARIES})
endif()


# Set common libs. Target specific libs set in subdirectories.
list(APPEND RMGLIBS RmgLib)
list(APPEND RMGLIBS Misc)
list(APPEND RMGLIBS XC)
list(APPEND RMGLIBS dftd3)
list(APPEND RMGLIBS MG)
list(APPEND RMGLIBS Force)
list(APPEND RMGLIBS Input)
list(APPEND RMGLIBS Finite_diff)
list(APPEND RMGLIBS US_PP)
list(APPEND RMGLIBS zfp)
list(APPEND RMGLIBS symspg_static)
list(APPEND RMGLIBS Diagonalization)



# Common include files
include_directories(${Boost_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/Headers/")
include_directories("${PROJECT_BINARY_DIR}/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/RmgLib/")
include_directories("${PROJECT_SOURCE_DIR}/RmgLib/include")
include_directories("${PROJECT_SOURCE_DIR}/zfp/include")
include_directories("${PROJECT_SOURCE_DIR}/spglib/src")


if(QMCPACK_SUPPORT)
    include_directories("${PROJECT_SOURCE_DIR}/Interfaces/")
endif()


# Common code subdirectories
add_subdirectory (RmgLib EXCLUDE_FROM_ALL) 
add_subdirectory (Finite_diff EXCLUDE_FROM_ALL) 
add_subdirectory (Force EXCLUDE_FROM_ALL) 
add_subdirectory (Input EXCLUDE_FROM_ALL) 
add_subdirectory (MG EXCLUDE_FROM_ALL) 
add_subdirectory (Misc EXCLUDE_FROM_ALL) 
add_subdirectory (US_PP EXCLUDE_FROM_ALL) 
add_subdirectory (zfp EXCLUDE_FROM_ALL) 
add_subdirectory (spglib EXCLUDE_FROM_ALL) 
add_subdirectory (XC EXCLUDE_FROM_ALL) 
add_subdirectory (dftd3 EXCLUDE_FROM_ALL) 
add_subdirectory (InternalPseudo EXCLUDE_FROM_ALL) 
add_subdirectory (Diagonalization EXCLUDE_FROM_ALL) 
if(RMG_CUDA_ENABLED)
    add_subdirectory(Gpufuncs EXCLUDE_FROM_ALL)
endif()

if(RMG_HIP_ENABLED)
    add_subdirectory(Gpufuncs EXCLUDE_FROM_ALL)
endif()

# Target specific code subdirectories. If no target
# is given the base code rmg-cpu or rmg-gpu is built by default
add_subdirectory (RMG) 
add_subdirectory (RMG/Common EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Subdiag EXCLUDE_FROM_ALL)

add_subdirectory (ON EXCLUDE_FROM_ALL) 
add_subdirectory (ON/Common EXCLUDE_FROM_ALL) 
add_subdirectory (ON/Input EXCLUDE_FROM_ALL) 
add_subdirectory (ON/ON-NEGF-share EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/ELDYN EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/IO_TDDFT EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/RMG_TDDFT EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/ON_TDDFT EXCLUDE_FROM_ALL) 

add_subdirectory (NEGF EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Common EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/SelfEnergy EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Input EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Local_current EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Pois_for_negf EXCLUDE_FROM_ALL) 

if(QMCPACK_SUPPORT)
    add_subdirectory (Interfaces EXCLUDE_FROM_ALL)
endif()
if(NOT CRAY_HOST)
    if(USE_INTERNAL_SCALAPACK)
        add_subdirectory (scalapack EXCLUDE_FROM_ALL) 
    else()
        if(NOT SCALAPACK_FOUND)
            add_subdirectory (scalapack EXCLUDE_FROM_ALL) 
	endif()
    endif()
endif()



#install(TARGETS ${rmg} 
#DESTINATION bin)

# This is mainly for packages since CPack will pick these up and install them
install(DIRECTORY Examples DESTINATION share/rmg/)
install(DIRECTORY tests DESTINATION share/rmg/)


# Packaging section next
set(cpack_input_file "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
set(CPACK_INPUT_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
set(CPACK_PROJECT_CONFIG_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
set(CPACK_GENERATOR "RPM")

set(CPACK_PACKAGE_VERSION_MAJOR "4")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(COMPLETE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})
set(CPACK_PACKAGE_NAME "rmg")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Emil Briggs")
set(CPACK_PACKAGE_VENDOR "North Carolina State University")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

include(CPack)

include(${CMAKE_SOURCE_DIR}/CTestConfig.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/macro.cmake)

enable_testing()
include(CTest)
SUBDIRS(tests)
