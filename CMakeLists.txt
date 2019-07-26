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
#        rmg-tddft-cpu  TDDFT based on rmg-on cpu only
#        rmg-tddft-gpu  TDDFT based on rmg-on with gpu support

#
# Environment variables that may also control configuration
#
# MKLROOT = root directory of Intel MKL libraries if installed.
# BOOST_ROOT = root directory of Boost installation
# CC  = Full path to C compiler
# CXX = Full path to C++ compiler
# FC  = Full path to Fortran compiler
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

include(FindSubversion)
if(Subversion_FOUND AND USE_SVN)
    execute_process (
        COMMAND svn info .. OUTPUT_VARIABLE foo)
    if(foo)
        string(REPLACE "Last Changed" ";" foo2 ${foo})
        list(GET foo2 2 foo3)
        string(REPLACE "\n" ";" foo4 ${foo3})
        string(REPLACE " " ";" foo3 ${foo4})
        list(GET foo3 2 CodeRevision)

        list(GET foo2 3 foo3)
        string(REPLACE " " ";" foo4 ${foo3})
        list(GET foo4 2 LastChangeDate)
    elseif()
        set(CodeRevision "release 3.0")
        set(LastChangeDate "Unknown")
    endif()
else(Subversion_FOUND AND USE_SVN)
    set(CodeRevision "3.0")
    set(LastChangeDate "Unknown")
endif(Subversion_FOUND AND USE_SVN)
add_definitions(-DSVN_REV="${CodeRevision}")
add_definitions(-DSVN_REVDATE="${LastChangeDate}")


# Require a 64 bit system
if (NOT CMAKE_SIZEOF_VOID_P MATCHES 8 )
    message( SEND_ERROR "RMG requires a 64 bit architechture." )
endif (NOT CMAKE_SIZEOF_VOID_P MATCHES 8)


# Check if this is a cray
if(DEFINED ENV{CRAY_LIBSCI_VERSION})
    # Cray platforms
    set (CRAY_XK 1)
endif()


# User settable options from cmake-gui, ccmake for text terminals
# or command line
option(RMG_GPU_ENABLED "Enable GPU support." OFF)
option(USE_MAGMA_LIBS "Use magma libraries. Requires RMG_GPU_ENABLED." OFF)
option(USE_OPENBABEL_LIBS "Openbabel support for importing atomic data." OFF)
option(USE_PLPLOT_LIBS "Use plplot for runtime graphs." OFF)
option(USE_INTERNAL_PSEUDOPOTENTIALS "Include internal pseudopotentials in binary." ON)
option(QMCPACK_SUPPORT "Generates a QMCPACK restart file." OFF)

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
option(USE_LIBXC "Use external libxc. Not recommended in most cases." OFF)
option(USE_SVN "Use svn to include version info in binary." ON)

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)



# Check for c++14 and c99 support
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    set(CMAKE_CXX_STANDARD 11)
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set(CMAKE_CXX_STANDARD 11)
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

# Send the Fortran C interface mangling template file to the build directory and set it up there
configure_file (
  "${PROJECT_SOURCE_DIR}/rmg_mangling.h.in"
  "${PROJECT_BINARY_DIR}/Headers/rmg_mangling.h"
)
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/Headers/rmg_mangling.h MACRO_NAMESPACE RMG_FC_)

# Here set any extra compiler flags required for a specific platform
# This is for stampede
#set(RMG_EXTRA_COMPILER_FLAGS "-gcc-name=/opt/apps/gcc/4.9.1/bin/gcc -xHost")

# Ryzen Threadripper with gcc
#
#set(RMG_EXTRA_COMPILER_FLAGS "-march=znver1 -Ofast -funroll-loops")
#set(RMG_EXTRA_COMPILER_FLAGS "-funroll-loops")
# Generic debug
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)


# Load custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Search for an MPI implementation. On Crays cmake does not pick up some
# vars but the compiler wrappers take care of it automatically.
if(CRAY_XK)
    find_package(MPI)
    CMAKE_FORCE_CXX_COMPILER(CC ${CMAKE_CXX_COMPILER_ID})
    CMAKE_FORCE_C_COMPILER(cc ${CMAKE_C_COMPILER_ID})
else()
    find_package(MPI REQUIRED)
endif()

# Get Openmp flags
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

if(USE_INTERNAL_PSEUDOPOTENTIALS)
    add_definitions(-DINTERNAL_PP=1)
endif()


if(RMG_GPU_ENABLED)
    add_definitions(-DGPU_ENABLED=1)
    set(rmg-negf "rmg-negf-gpu")
    set(rmg "rmg-gpu")
    set(rmg-on "rmg-on-gpu")
    set(rmg-tddft "rmg-tddft-gpu")
else(RMG_GPU_ENABLED)
    add_definitions(-DGPU_ENABLED=0)
    add_definitions(-DMAGMA_LIBS=0)
    set(rmg-negf "rmg-negf-cpu")
    set(rmg "rmg-cpu")
    set(rmg-on "rmg-on-cpu")
    set(rmg-tddft "rmg-tddft-cpu")
endif(RMG_GPU_ENABLED)


# System specific flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(OperatingSystem "Linux")
    add_definitions(-DLINUX=1)
    set(ZLIB "z")
    set(RTLIB "rt")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")


# GNUCC specific compiler flags
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pthread -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -Wno-write-strings -pthread -m64 -fext-numeric-literals")
    set (CMAKE_Fortran_FLAGS_RELEASE "-lg2c -cpp -c")
    set (CMAKE_Fortran_FLAGS "-lg2c -cpp -c")
    set(CMAKE_EXE_LINKER_FLAGS "-ldl")
    set(RMG_EXTRA_COMPILER_FLAGS "${RMG_EXTRA_COMPILER_FLAGS} -fcx-limited-range -fno-trapping-math")
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

if(USE_LIBXC)
    add_definitions(-D__LIBXC=1)
endif()

if (CRAY_XK)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lhugetlbfs")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lhugetlbfs")
   # External packages sometimes have trouble with this on Crays
   set($ENV{MPICC} "${CMAKE_C_COMPILER}")
endif()


# This file contains system specific paths. If cmake is unable to find
# a specific package or library look in this file and uncomment or add
# specific search paths.
include(CMakeFindRootPath.inc)


# Executables defined by the file that contains main()
# Main codes
add_executable(${rmg} RMG/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
add_executable(${rmg-tddft} TDDFT/Common/Main_rmg.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
add_executable(${rmg-on} ON/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
add_executable(on-tddft TDDFT/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
add_executable(${rmg-negf} NEGF/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)

# adds qmcpack support to base code
if(QMCPACK_SUPPORT)
    add_definitions(-DQMCPACK_SUPPORT=1)
    target_link_libraries(${rmg} Interfaces RMG_Common Interfaces)
endif()


if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    target_link_libraries (${rmg} ${MPI_CXX_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${MPI_CXX_LIBRARIES})
    target_link_libraries (${rmg-on} ${MPI_CXX_LIBRARIES})
    target_link_libraries (on-tddft ${MPI_CXX_LIBRARIES})
    target_link_libraries (${rmg-negf} ${MPI_CXX_LIBRARIES})
endif()

if(MPI_C_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
    target_link_libraries (${rmg} ${MPI_C_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${MPI_C_LIBRARIES})
    target_link_libraries (${rmg-on} ${MPI_C_LIBRARIES})
    target_link_libraries (on-tddft ${MPI_C_LIBRARIES})
    target_link_libraries (${rmg-negf} ${MPI_C_LIBRARIES})
endif()


# Require pthreads on UNIX
if (${UNIX})
    find_package(Threads REQUIRED)
    target_link_libraries (${rmg} ${CMAKE_THREAD_LIBS_INIT} ${RTLIB})
    target_link_libraries (${rmg-tddft} ${CMAKE_THREAD_LIBS_INIT} ${RTLIB})
    target_link_libraries (${rmg-on} ${RTLIB})
endif(${UNIX})


find_package(NUMA)
if(NUMA_FOUND AND NUMA_INCLUDES)
    include_directories(${NUMA_INCLUDES})
    add_definitions(-DUSE_NUMA=1)
    set(USE_NUMA_INTERNAL 1)
else(NUMA_FOUND AND NUMA_INCLUDES)
    set(USE_NUMA_INTERNAL 0)
endif(NUMA_FOUND AND NUMA_INCLUDES)


find_package(HWLOC)
if(HWLOC_FOUND AND HWLOC_INCLUDES)
    include_directories(${HWLOC_INCLUDES})
    #add_definitions(-DUSE_HWLOC=1)
    target_link_libraries (${rmg} ${HWLOC_LIBRARIES})
    target_link_libraries (${rmg-on} ${HWLOC_LIBRARIES})
    target_link_libraries (on-tddft ${HWLOC_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${HWLOC_LIBRARIES})
    target_link_libraries (${rmg-negf} ${HWLOC_LIBRARIES})
endif(HWLOC_FOUND AND HWLOC_INCLUDES)


if(USE_ELPA_LIBS)
    add_definitions(-DUSE_ELPA=1)
    target_link_libraries(${rmg} /home/briggs/src/elpa-2016.11.001.pre/.libs/libelpa.so)
    #target_link_libraries(${rmg} /usr/lib64/libelpa.so)
endif(USE_ELPA_LIBS)


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "^(FreeBSD|DragonFly)$")
    if(CRAY_XK)
        add_definitions(-DHAVE_ASYNC_ALLREDUCE=1)
    else(CRAY_XK)
        include("PlatformChecks/CheckAsyncAllreduce.inc")
        if(RMG_ASYNC_ALLREDUCE_INTERNAL)
            add_definitions(-DHAVE_ASYNC_ALLREDUCE=1)
        else(RMG_ASYNC_ALLREDUCE_INTERNAL)
            add_definitions(-DHAVE_ASYNC_ALLREDUCE=0)
        endif(RMG_ASYNC_ALLREDUCE_INTERNAL)
    endif(CRAY_XK)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "^(FreeBSD|DragonFly)$")


if(USE_PLPLOT_LIBS)
    find_package(PLPLOT)
    if(PLPLOT_FOUND)
        include("PlatformChecks/CheckPlplot.inc")
        if(USE_PLPLOT_LIBS_INTERNAL)
            include_directories(${PLplot_INCLUDE_DIR})
            include_directories("${PLplot_INCLUDE_DIR}/plplot")
            add_definitions(-DPLPLOT_LIBS=1)
        else(USE_PLPLOT_LIBS_INTERNAL)
            add_definitions(-DPLPLOT_LIBS=0)
        endif(USE_PLPLOT_LIBS_INTERNAL)
    else(PLPLOT_FOUND)
        add_definitions(-DPLPLOT_LIBS=0)
    endif(PLPLOT_FOUND)
endif(USE_PLPLOT_LIBS)


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(BOOST_ROOT "...")
    if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT  $ENV{BOOST_ROOT})
    endif()
    find_package(Boost REQUIRED COMPONENTS thread system iostreams program_options filesystem REQUIRED)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")


if(CRAY_XK)
    find_package(BZ2 REQUIRED)
    find_package(Z REQUIRED)
    target_link_libraries (${rmg} ${BZ2_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${BZ2_LIBRARIES})
    target_link_libraries (${rmg-on} ${BZ2_LIBRARIES})
    target_link_libraries (on-tddft ${BZ2_LIBRARIES})
    target_link_libraries (${rmg-negf} ${BZ2_LIBRARIES})
    target_link_libraries (${rmg} ${Z_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${Z_LIBRARIES})
    target_link_libraries (${rmg-on} ${Z_LIBRARIES})
    target_link_libraries (on-tddft ${Z_LIBRARIES})
    target_link_libraries (${rmg-negf} ${Z_LIBRARIES})
endif(CRAY_XK)


# Look for openbabel
if(USE_OPENBABEL_LIBS)
    find_package(OPENBABEL)
    include("PlatformChecks/CheckOpenbabel.inc")

    if(CHECK_OPENBABEL_INTERNAL)
        include_directories(${OPENBABEL_INCLUDES})
        add_definitions(-DOPENBABEL_LIBS=1)
    else(CHECK_OPENBABEL_INTERNAL)
        add_definitions(-DOPENBABEL_LIBS=0)
    endif (CHECK_OPENBABEL_INTERNAL)

endif(USE_OPENBABEL_LIBS)


# Look for the magma libraries
if(USE_MAGMA_LIBS AND RMG_GPU_ENABLED)
  find_package(MAGMA)
  if(MAGMA_FOUND)
    add_definitions(-DMAGMA_LIBS=1)
    include_directories(${MAGMA_INCLUDES})
    if(MAGMABLAS_FOUND) 
        target_link_libraries (${rmg} ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (${rmg-tddft} ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (${rmg-on} ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (on-tddft ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (${rmg-negf} ${MAGMA_LIBRARIES} ${MAGMABLAS_LIBRARIES} )
    else(MAGMABLAS_FOUND) 
        target_link_libraries (${rmg} ${MAGMA_LIBRARIES} )
        target_link_libraries (${rmg-tddft} ${MAGMA_LIBRARIES} )
        target_link_libraries (${rmg-on} ${MAGMA_LIBRARIES} )
        target_link_libraries (on-tddft ${MAGMA_LIBRARIES} )
        target_link_libraries (${rmg-negf} ${MAGMA_LIBRARIES} )
    endif(MAGMABLAS_FOUND) 
  elseif(MAGMA_FOUND)
      add_definitions(-DMAGMA_LIBS=0)
  endif(MAGMA_FOUND)
endif(USE_MAGMA_LIBS AND RMG_GPU_ENABLED)


# For various math libraries start by looking for Intel MKL. No good cmake module currently for MKL.
# Check if environment variable is set
if(NOT MKLROOT)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # Check the default location
        if (EXISTS "/opt/intel/mkl")
#            set(MKLROOT "/opt/intel/mkl")
        endif (EXISTS "/opt/intel/mkl")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
endif()


if(NOT CRAY_XK)
    if(MKLROOT)
        include_directories(${MKLROOT}/include)
        add_definitions(-DUSE_MKL)
        set(FFTW_INCLUDES "${MKLROOT}/include/fftw")
        set(MKL_LIBDIR "${MKLROOT}/lib/intel64")
        set(MKL_LIB64 "${MKLROOT}/lib/intel64")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            if(EXISTS "${MKL_LIB64}")
                set(USE_LINUX_MKL_LIBS 1)
                target_link_libraries (${rmg} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
                target_link_libraries (${rmg-tddft} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
                target_link_libraries (${rmg-on} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
                target_link_libraries (${rmg-negf} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
                set(RMG_MATH_LIBS "done")
            else(EXISTS "${MKL_LIB64}")
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif(EXISTS "${MKL_LIB64}")
        endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    endif(MKLROOT)
endif(NOT CRAY_XK)


# MKL includes fftw blas, blacs, lapack and scalapack. If it is not present we have to look for these libraries
# individually
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
target_link_libraries (${rmg} ${FFTW_LIBRARIES})
target_link_libraries (${rmg-tddft} ${FFTW_LIBRARIES})
target_link_libraries (${rmg-on} ${FFTW_LIBRARIES})
target_link_libraries (on-tddft ${FFTW_LIBRARIES})
target_link_libraries (${rmg-negf} ${FFTW_LIBRARIES})
target_link_libraries (${rmg} RMG_Common ${FFTW_LIBRARIES})
target_link_libraries (${rmg-tddft} RMG_Common ${FFTW_LIBRARIES})


if(NOT RMG_MATH_LIBS)

    # FFTW
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDES})
    target_link_libraries (${rmg} ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-on} ${FFTW_LIBRARIES})
    target_link_libraries (on-tddft ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-negf} ${FFTW_LIBRARIES})
    target_link_libraries (${rmg} RMG_Common ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-tddft} RMG_Common ${FFTW_LIBRARIES})

    if(NOT CRAY_XK)

        # BLAS
        find_package(BLAS REQUIRED)
        if(BLAS_FOUND)
            target_link_libraries (${rmg} ${BLAS_LIBRARIES})
            target_link_libraries (${rmg-tddft} ${BLAS_LIBRARIES})
            target_link_libraries (${rmg-on} ${BLAS_LIBRARIES})
            target_link_libraries (on-tddft ${BLAS_LIBRARIES})
            target_link_libraries (${rmg-negf} ${BLAS_LIBRARIES})
        endif (BLAS_FOUND)

    endif (NOT CRAY_XK)

endif(NOT RMG_MATH_LIBS)


# Check for openblas library supporting changing threads
include("PlatformChecks/CheckOpenBlasSetNumThreads.inc")
if(RMG_OPENBLAS_SET_NUM_THREADS_INTERNAL)
    add_definitions(-DOPENBLAS_SET_NUM_THREADS=1)
else(RMG_OPENBLAS_SET_NUM_THREADS_INTERNAL)
    add_definitions(-DOPENBLAS_SET_NUM_THREADS=0)
endif(RMG_OPENBLAS_SET_NUM_THREADS_INTERNAL)


# HDF5 libraries for QMCPACK integration
SET(ENABLE_PHDF5 1 CACHE BOOL "Enable phdf5 for output")
SET(HAVE_LIBHDF5 0)
SET(HDF5_PREFER_PARALLEL 0)
IF(HAVE_MPI)
 SET(HDF5_PREFER_PARALLEL 1)
ENDIF()
find_package(HDF5 COMPONENTS C)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  IF ( CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG )
    IF ( CMAKE_BUILD_TYPE MATCHES DEBUG )
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG} )
    ELSE()
      SET( HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE} )
    ENDIF()
  ENDIF()
  IF(HDF5_IS_PARALLEL)
    IF(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 parallel collective I/O")
    ELSE(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O")
    ENDIF(ENABLE_PHDF5)
  ELSE(HDF5_IS_PARALLEL)
    IF(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O due to the lack of library support for parallel")
      SET(ENABLE_PHDF5 0)
    ELSE(ENABLE_PHDF5)
      MESSAGE(STATUS "Using HDF5 non-scalable serial I/O")
    ENDIF(ENABLE_PHDF5)
  ENDIF(HDF5_IS_PARALLEL)
ELSE(HDF5_FOUND)
#    MESSAGE(FATAL_ERROR "Require hdf5 1.6.4 or higher. Set HDF5_ROOT")
ENDIF(HDF5_FOUND)



# Look for GPU libraries. Only CUDA supported for now
if(RMG_GPU_ENABLED)
  find_package(CUDA REQUIRED)
  if(${CUDA_VERSION_MAJOR} LESS 6)
      message(FATAL_ERROR "\nCompiling RMG with GPU support requires CUDA version 6.0 or greater. Either set GPU_ENABLED=0 to build a CPU only version or upgrade your cuda installation.") 
  endif(${CUDA_VERSION_MAJOR} LESS 6)

  if(NOT CRAY_XK)
      target_link_libraries(${rmg} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY})
      target_link_libraries(${rmg-tddft} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY})
      target_link_libraries(${rmg-on} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY})
      target_link_libraries(on-tddft ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY})
      target_link_libraries(${rmg-negf} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
  endif(NOT CRAY_XK)
  target_link_libraries(${rmg} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(${rmg-tddft} ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(${rmg-on} ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY})
  target_link_libraries(${rmg-on} ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY})
  target_link_libraries(on-tddft ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(on-tddft ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(${rmg-negf} ${CUDA_CUDART_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(${rmg-negf} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
  include_directories(${CUDA_INCLUDE_DIRS})


# Cuda functions that we build are in Gpufuncs/CMakeLists.txt
 
  add_definitions(-DCUDA_VERSION_MAJOR=${CUDA_VERSION_MAJOR})
  add_definitions(-DCUDA_VERSION_MINOR=${CUDA_VERSION_MINOR})
endif(RMG_GPU_ENABLED)

#include_directories(rmg "/usr/local/magma/include")

if(USE_LIBXC)
    find_package(LIBXC)
    if(NOT LIBXC_FOUND)
        execute_process(
        COMMAND tar xvf ${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2.tar.gz -C ${PROJECT_SOURCE_DIR}/lib )
        execute_process(
        COMMAND mkdir ${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64)
        execute_process(
        COMMAND ${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/configure
            --prefix=${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64
            --libdir=${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64)
        execute_process(COMMAND make install
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64)
    endif(NOT LIBXC_FOUND)
    find_package(LIBXC REQUIRED)
    include_directories("${LIBXC_INCLUDES}"  "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/")
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
    endif(NOT ELEMENTAL_FOUND)

    find_package(ELEMENTAL REQUIRED)
    include_directories(${ELEMENTAL_INCLUDES})
    target_link_libraries (${rmg-on} ON_Common ${ELEMENTAL_LIBRARIES})
    find_package(PMRRR)
    include_directories(${PMRRR_INCLUDES})
    target_link_libraries (${rmg-on} ${PMRRR_LIBRARIES})
    target_link_libraries (${rmg-on} ON_Common ${PMRRR_LIBRARIES})
endif()


# include files
include_directories(${Boost_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/Headers/")
include_directories("${PROJECT_BINARY_DIR}/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/RmgLib/")
include_directories("${PROJECT_SOURCE_DIR}/RmgLib/include")
include_directories("${PROJECT_SOURCE_DIR}/zfp/include")
include_directories("${PROJECT_SOURCE_DIR}/spglib/src")
include_directories("${PROJECT_SOURCE_DIR}/ON/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/TDDFT/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/NEGF/Headers/")
if(QMCPACK_SUPPORT)
    include_directories("${PROJECT_SOURCE_DIR}/Interfaces/")
endif()


# Common code subdirectories
add_subdirectory (RmgLib) 
add_subdirectory (Finite_diff) 
add_subdirectory (Force) 
add_subdirectory (Input) 
add_subdirectory (MG) 
add_subdirectory (Misc) 
add_subdirectory (US_PP) 
add_subdirectory (zfp) 
add_subdirectory (spglib) 
if(QMCPACK_SUPPORT)
    add_subdirectory (Interfaces)
endif()
if(NOT CRAY_XK)
    add_subdirectory (scalapack) 
endif()
if(USE_LIBXC)
    add_subdirectory (XC_useLIBXC) 
endif(USE_LIBXC)
add_subdirectory (XC) 
add_subdirectory (InternalPseudo) 
if(RMG_GPU_ENABLED)
    add_subdirectory(Gpufuncs)
endif(RMG_GPU_ENABLED)

# Target specific code subdirectories
add_subdirectory (RMG/Common EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Subdiag EXCLUDE_FROM_ALL)

add_subdirectory (ON/Common EXCLUDE_FROM_ALL) 
add_subdirectory (ON/Input EXCLUDE_FROM_ALL) 
add_subdirectory (ON/ON-NEGF-share EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/Common EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/ELDYN EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/OrbitalProjection EXCLUDE_FROM_ALL) 

add_subdirectory (NEGF/Common EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/SelfEnergy EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Input EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Local_current EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Pois_for_negf EXCLUDE_FROM_ALL) 

# Link subdirectories
link_directories (${RMG_BINARY_DIR}/Finite_diff) 
link_directories (${RMG_BINARY_DIR}/Force)
link_directories (${RMG_BINARY_DIR}/Input)
link_directories (${RMG_BINARY_DIR}/MG)
link_directories (${RMG_BINARY_DIR}/Misc)
link_directories (${RMG_BINARY_DIR}/zfp)
link_directories (${RMG_BINARY_DIR}/spglib)
if(QMCPACK_SUPPORT)
    link_directories (${RMG_BINARY_DIR}/Interfaces)
endif()
if(NOT CRAY_XK)
    link_directories (${RMG_BINARY_DIR}/scalapack)
endif()
link_directories (${RMG_BINARY_DIR}/US_PP)
link_directories (${RMG_BINARY_DIR}/XC)

link_directories (${RMG_BINARY_DIR}/RMG/Common)
link_directories (${RMG_BINARY_DIR}/RMG/Input)
link_directories (${RMG_BINARY_DIR}/RMG/Spin)
link_directories (${RMG_BINARY_DIR}/RMG/Spin/XC)
link_directories (${RMG_BINARY_DIR}/RMG/Subdiag)


# These include files include some subproject specific stuff
# in order to keep from cluttering up the main CMakeLists.txt
# file.
include("SubprojectIncludes/CMakeLists.txt.rmg.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-tddft.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-on.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-negf.inc")
include("SubprojectIncludes/CMakeLists.txt.on-tddft.inc")
if(USE_LINUX_MKL_LIBS)
    target_link_libraries (${rmg-tddft} -Wl,--start-group ${MKL_LIB64}/libfftw3x_cdft_lp64.a ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
    target_link_libraries (${rmg} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group)
    target_link_libraries (${rmg-on} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group)
    target_link_libraries (${rmg-negf} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group)
endif(USE_LINUX_MKL_LIBS)

install(TARGETS ${rmg} 
DESTINATION bin)

#install(FILES RmgLib/libRmgLibShared.so RmgLib/libRmgLib.a
#DESTINATION lib)

install(DIRECTORY Examples DESTINATION share/rmg)

# Packaging section next
SET(cpack_input_file "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
SET(CPACK_INPUT_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
SET(CPACK_PROJECT_CONFIG_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
set(CPACK_GENERATOR "RPM")

set(CPACK_PACKAGE_VERSION_MAJOR "4")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(COMPLETE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})
set(CPACK_PACKAGE_NAME "rmg")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Emil Briggs")
set(CPACK_PACKAGE_VENDOR "North Carolina State University")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

include(CPack)

