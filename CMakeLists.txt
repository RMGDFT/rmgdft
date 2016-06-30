#
#    RMG master cmake file
#
#    Subprojects:
#        rmg        Base code
#        rmg-on     ON code
#        NEGF
#        on-tddft:  TDDFT based on rmg-on
#
# Environment variables that may also control configuration
# MKLROOT = root directory of Intel MKL libraries if installed.
# BOOST_ROOT = root directory of Boost installation
# CC = Full path to C compiler
# CXX = Full path to C++ compiler
#
#
cmake_minimum_required (VERSION 2.8) 

enable_language(Fortran)
project (RMG C CXX Fortran) 

# Standard modules for platform checks
INCLUDE (CheckFunctionExists)
INCLUDE (CheckSymbolExists)
INCLUDE (CheckCSourceCompiles)
INCLUDE (CheckCXXSourceCompiles)
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckLibraryExists)
INCLUDE (CMakeForceCompiler)


# User editable options
set (RMG_GPU_ENABLED 1)
set (RMG_PFFT_LIBS 1) 
set (RMG_MAGMA_LIBS 1) 
set (USE_SCALAPACK_LIBS 1)
set (USE_ELEMENTAL_LIBS 0)
set (USE_LIBXC 0)
set (USE_OPENBABEL_LIBS 1)
set (USE_SVN 0)

# Load custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)


if(RMG_GPU_ENABLED)
    set(rmg-negf "rmg-negf-gpu")
    set(rmg "rmg-gpu")
    set(rmg-tddft "rmg-tddft-gpu")
else(RMG_GPU_ENABLED)
    set(rmg-negf "rmg-negf-cpu")
    set(rmg "rmg-cpu")
    set(rmg-tddft "rmg-tddft-cpu")
endif(RMG_GPU_ENABLED)

# Here set any extra compiler flags required for a specific platform
# This is for stampede
#set(RMG_EXTRA_COMPILER_FLAGS "-gcc-name=/opt/apps/gcc/4.9.1/bin/gcc -xHost")


# Detect system type and set appropriate flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    set(OperatingSystem "Linux")
    set(RMG_SYSTEM_COMPILER_FLAG "-DLINUX=1")
    set(ZLIB "z")
    set(RTLIB "rt")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Windows specific code. RMG currently only builds on Windows with
    # the Intel compilers and MKL math library. Prebuilt versions of
    # the libxc and spglib libraries are also supplied in the source distribtion
    set(OperatingSystem "Windows")
    set(MKLROOT "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/")
    set(CMAKE_Fortran_COMPILER "C:/Program Files (x86)/Intel/Composer XE 2015/bin/intel64/ifort.exe")
    set(CMAKE_C_COMPILER "C:/Program Files (x86)/Intel/Composer XE 2015/bin/intel64/icl.exe")
    set(CMAKE_CXX_COMPILER "C:/Program Files (x86)/Intel/Composer XE 2015/bin/intel64/icl.exe")
    set(CMAKE_CXX_FLAGS "/DWIN32 /D_WIN32 /D_WINDOWS /W3 /GR /EHsc /Qstd=c++11 /DNOMINMAX /D_NOMINMAX /Qopenmp" CACHE STRING "" FORCE)
    set(CMAKE_C_FLAGS "/DWIN32 /D_WIN32 /D_WINDOWS /W3 /Qstd=c99 /DNOMINMAX /D_NOMINMAX /Qopenmp" CACHE STRING "" FORCE)

    set(LIBXC_LIBRARIES "${PROJECT_SOURCE_DIR}/lib/winlibs/x64/LIBXC.lib")
    set(LIBSPG_LIBRARIES "${PROJECT_SOURCE_DIR}/lib/winlibs/x64/SPGLIBVC12.lib")    
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(DEFINED ENV{CRAY_LIBSCI_VERSION})
    # Cray platforms
    set (CRAY_XK 1)
endif()


# Compiler flags. We currently only support GNU and Intel compilers so if
# GNUCC is not detected we set the compiler to Intel
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "-Wall -O3 ${RMG_SYSTEM_COMPILER_FLAG} ${RMG_EXTRA_COMPILER_FLAGS} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -mfpmath=sse -fopenmp -m64 -I$(MKLROOT)/include")
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs -O3 -std=c++11 -Wno-write-strings ${RMG_SYSTEM_COMPILER_FLAG} ${RMG_EXTRA_COMPILER_FLAGS} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -mfpmath=sse -fopenmp -m64 -pthread -I$(MKLROOT)/include")
    set(CMAKE_EXE_LINKER_FLAGS "-ldl -fopenmp")
    if(USE_LIBXC)
        set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c -cpp -D__LIBXC=1 -c")
        set (CMAKE_Fortran_FLAGS "-O3 -fopenmp -march=native -lg2c -cpp -D__LIBXC=1 -c")
    else(USE_LIBXC)
        set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c -cpp -c")
        set (CMAKE_Fortran_FLAGS "-O3 -fopenmp -march=native -lg2c -cpp -c")
    endif(USE_LIBXC)
    set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -static-libstdc++ -ldl")
else(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # Intel compiler
    set(CMAKE_CXX_FLAGS "${RMG_EXTRA_COMPILER_FLAGS} -O3 -std=c++0x -DBATCH_NLS=1 -DNOMINMAX -D_NOMINMAX -openmp -I${MKLROOT}/include ${RMG_SYSTEM_COMPILER_FLAG}")
    set(CMAKE_C_FLAGS "${RMG_EXTRA_COMPILER_FLAGS} -O3 -std=c99 -DBATCH_NLS=1 -DNOMINMAX -D_NOMINMAX -openmp -I${MKLROOT}/include ${RMG_SYSTEM_COMPILER_FLAG}")
    set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -static-libstdc++")
endif()


if (CRAY_XK)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lhugetlbfs")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lhugetlbfs")
   set(USE_SCALAPACK_LIBS 1)
   # External packages (e.g. PFFT) need this
   set($ENV{MPICC} "${CMAKE_C_COMPILER}")
endif(CRAY_XK)


if (RMG_GPU_ENABLED)
   set (GPU_ENABLED 1)
   if(RMG_MAGMA_LIBS)
      set (MAGMA_LIBS 1)
   else (RMG_MAGMA_LIBS)
      set (MAGMA_LIBS 0)
   endif(RMG_MAGMA_LIBS) 
else (RMG_GPU_ENABLED)
   set (GPU_ENABLED 0)
   set (MAGMA_LIBS 0)
endif(RMG_GPU_ENABLED)

include(CMakeFindRootPath.inc)
# Executables defined by the file that contains main()
# Main codes
add_executable(${rmg} RMG/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(${rmg-tddft} TDDFT/Common/Main_rmg.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-on ON/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# ON code
add_executable(on-tddft TDDFT/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)
add_executable(${rmg-negf} NEGF/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)


# Require a 64 bit system
if (NOT CMAKE_SIZEOF_VOID_P MATCHES 8 )
    message( SEND_ERROR "RMG requires a 64 bit architechture." )
endif (NOT CMAKE_SIZEOF_VOID_P MATCHES 8)



# Search for an MPI implementation
#set(MPI_C_COMPILER "/opt/openmpi/bin/mpicc")
#find_package(MPI REQUIRED)
find_package(MPI)
#  include_directories(${MPI_INCLUDE_PATH})
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    target_link_libraries (${rmg} ${MPI_CXX_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${MPI_CXX_LIBRARIES})
    target_link_libraries (rmg-on ${MPI_CXX_LIBRARIES})
    target_link_libraries (on-tddft ${MPI_CXX_LIBRARIES})
    target_link_libraries (${rmg-negf} ${MPI_CXX_LIBRARIES})
else (MPI_CXX_FOUND)
    if(CRAY_XK)
	CMAKE_FORCE_CXX_COMPILER(CC ${CMAKE_CXX_COMPILER_ID})
    endif(CRAY_XK)	
endif (MPI_CXX_FOUND)
if(MPI_C_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
    target_link_libraries (${rmg} ${MPI_C_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${MPI_C_LIBRARIES})
    target_link_libraries (rmg-on ${MPI_C_LIBRARIES})
    target_link_libraries (on-tddft ${MPI_C_LIBRARIES})
    target_link_libraries (${rmg-negf} ${MPI_C_LIBRARIES})
else(MPI_C_FOUND)
    CMAKE_FORCE_C_COMPILER(cc ${CMAKE_C_COMPILER_ID})
endif (MPI_C_FOUND)


# Require pthreads on UNIX
if (${UNIX})
    set (CMAKE_THREAD_PREFER_PTHREAD)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
        target_link_libraries (${rmg} ${CMAKE_THREAD_LIBS_INIT} ${RTLIB})
        target_link_libraries (${rmg-tddft} ${CMAKE_THREAD_LIBS_INIT} ${RTLIB})
    endif(CMAKE_USE_PTHREADS_INIT)
    target_link_libraries (rmg-on ${RTLIB})
endif(${UNIX})



if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include("PlatformChecks/CheckAsyncAllreduce.inc")
    if(RMG_ASYNC_ALLREDUCE_INTERNAL)
        set(RMG_ASYNC_ALLREDUCE 1)
    else(RMG_ASYNC_ALLREDUCE_INTERNAL)
        set(RMG_ASYNC_ALLREDUCE 0)
    endif(RMG_ASYNC_ALLREDUCE_INTERNAL)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

find_package(PLPLOT)
if(PLPLOT_FOUND)
    include("PlatformChecks/CheckPlplot.inc")
    if(USE_PLPLOT_LIBS_INTERNAL)
        include_directories(${PLplot_INCLUDE_DIR})
        include_directories("${PLplot_INCLUDE_DIR}/plplot")
        set(USE_PLPLOT_LIBS 1)
    else(USE_PLPLOT_LIBS_INTERNAL)
        set(USE_PLPLOT_LIBS 0)
    endif(USE_PLPLOT_LIBS_INTERNAL)
else(PLPLOT_FOUND)
    set(USE_PLPLOT_LIBS 0)
endif(PLPLOT_FOUND)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(RMG_ASYNC_ALLREDUCE 0)
    #message("CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
    set(Boost_USE_STATIC_RUNTIME    OFF)
    set(Boost_COMPILER "-vc120")
    set(BOOST_ROOT "C:/local/boost_1_57_0/")
    set(Boost_DEBUG ON)
    set(BOOST_LIBRARYDIR "C:/local/boost_1_57_0/lib64-msvc-12.0/")
    set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0")
    find_package(Boost COMPONENTS thread system iostreams program_options filesystem date_time chrono zlib bzip2 REQUIRED)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(BOOST_ROOT "...")
    if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT  $ENV{BOOST_ROOT})
    endif()
    set(Boost_USE_STATIC_LIBS        OFF)
    set(Boost_USE_MULTITHREADED      ON)
    find_package(Boost REQUIRED COMPONENTS thread system iostreams program_options filesystem REQUIRED)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")



if(CRAY_XK)
    find_package(BZ2 REQUIRED)
    find_package(Z REQUIRED)
    target_link_libraries (${rmg} ${BZ2_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${BZ2_LIBRARIES})
    target_link_libraries (rmg-on ${BZ2_LIBRARIES})
    target_link_libraries (on-tddft ${BZ2_LIBRARIES})
    target_link_libraries (${rmg-negf} ${BZ2_LIBRARIES})
    target_link_libraries (${rmg} ${Z_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${Z_LIBRARIES})
    target_link_libraries (rmg-on ${Z_LIBRARIES})
    target_link_libraries (on-tddft ${Z_LIBRARIES})
    target_link_libraries (${rmg-negf} ${Z_LIBRARIES})
endif(CRAY_XK)


# Look for openbabel
find_package(OPENBABEL)
include("PlatformChecks/CheckOpenbabel.inc")

if(CHECK_OPENBABEL_INTERNAL)
    include_directories(${OPENBABEL_INCLUDES})
    unset (USE_OPENBABEL_LIBS CACHE)
    set(USE_OPENBABEL_LIBS 1)
else(CHECK_OPENBABEL_INTERNAL)
    unset (USE_OPENBABEL_LIBS CACHE)
    set(USE_OPENBABEL_LIBS 0)
endif (CHECK_OPENBABEL_INTERNAL)


# Look for the magma libraries
if(MAGMA_LIBS)
  find_package(MAGMA)
  include_directories(${MAGMA_INCLUDES})
  if(MAGMA_FOUND)
    if(MAGMABLAS_FOUND) 
        target_link_libraries (${rmg} ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (${rmg-tddft} ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (rmg-on ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (on-tddft ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (${rmg-negf} ${MAGMA_LIBRARIES} ${MAGMABLAS_LIBRARIES} )
    else(MAGMABLAS_FOUND) 
        target_link_libraries (${rmg} ${MAGMA_LIBRARIES} )
        target_link_libraries (${rmg-tddft} ${MAGMA_LIBRARIES} )
        target_link_libraries (rmg-on ${MAGMA_LIBRARIES} )
        target_link_libraries (on-tddft ${MAGMA_LIBRARIES} )
        target_link_libraries (${rmg-negf} ${MAGMA_LIBRARIES} )
    endif(MAGMABLAS_FOUND) 
  endif(MAGMA_FOUND)
endif(MAGMA_LIBS)



# For various math libraries start by looking for Intel MKL. No good cmake module currently for MKL.
# Check if environment variable is set
set(MKLROOT $ENV{MKLROOT})
if (NOT MKLROOT)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # Check the default location
        if (EXISTS "/opt/intel/mkl")
            set(MKLROOT "/opt/intel/mkl")
        endif (EXISTS "/opt/intel/mkl")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
endif (NOT MKLROOT)

if(NOT CRAY_XK)
    if(MKLROOT)
        set(FFTW_INCLUDES "${MKLROOT}/include/fftw")
        set(MKL_LIBDIR "${MKLROOT}/lib/intel64")
        set(MKL_LIB64 "${MKLROOT}/lib/intel64")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            if(EXISTS "${MKL_LIB64}")
                target_link_libraries (${rmg} "${MKL_LIB64}/mkl_intel_lp64.lib")
                target_link_libraries (${rmg-tddft} "${MKL_LIB64}/mkl_intel_lp64.lib")
                target_link_libraries (${rmg} "${MKL_LIB64}/mkl_core.lib")
                target_link_libraries (${rmg-tddft} "${MKL_LIB64}/mkl_core.lib")
                target_link_libraries (${rmg} "${MKL_LIB64}/mkl_intel_thread.lib")
                target_link_libraries (${rmg-tddft} "${MKL_LIB64}/mkl_intel_thread.lib")
                set(RMG_MATH_LIBS "done")
            else(EXISTS "${MKL_LIB64}")
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif(EXISTS "${MKL_LIB64}")
        else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            if(EXISTS "${MKL_LIB64}")
                set(USE_LINUX_MKL_LIBS 1)
                target_link_libraries (${rmg} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_scalapack_lp64.a ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a ${MKL_LIB64}/libmkl_blacs_openmpi_lp64.a -Wl,--end-group)
                target_link_libraries (${rmg-tddft} -Wl,--start-group ${MAGMA_LIBRARIES} ${MKL_LIB64}/libmkl_scalapack_lp64.a ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a ${MKL_LIB64}/libmkl_blacs_openmpi_lp64.a -Wl,--end-group)
                set(RMG_MATH_LIBS "done")
            else(EXISTS "${MKL_LIB64}")
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif(EXISTS "${MKL_LIB64}")
        endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    endif(MKLROOT)
endif(NOT CRAY_XK)

include_directories(${FFTW_INCLUDES})
# MKL includes fftw blas, blacs, lapack and scalapack. If it is not present we have to look for these libraries
# individually



if(NOT RMG_MATH_LIBS)

    # FFTW
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDES})
    target_link_libraries (${rmg} ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-tddft} ${FFTW_LIBRARIES})
    target_link_libraries (rmg-on ${FFTW_LIBRARIES})
    target_link_libraries (on-tddft ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-negf} ${FFTW_LIBRARIES})
    target_link_libraries (${rmg} RMG_Common ${FFTW_LIBRARIES})
    target_link_libraries (${rmg-tddft} RMG_Common ${FFTW_LIBRARIES})

    if(NOT CRAY_XK)

        # OPENBLAS
        set(OPEN_BLAS_SEARCH_PATHS  /lib/ /lib64/  /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /opt/OpenBLAS/lib $ENV{OPENBLAS_HOME}/lib )
        FIND_LIBRARY(OPENBLAS NAMES openblas PATHS ${OPEN_BLAS_SEARCH_PATHS})
        if (OPENBLAS)
            set(BLAS_FOUND TRUE)
            set(BLAS_LIBRARIES ${OPENBLAS})
            set(BLAS_blas_LIBRARY ${OPENBLAS})
        else (OPENBLAS)
            find_package(BLAS REQUIRED)
        endif (OPENBLAS)

        if(BLAS_FOUND)
            target_link_libraries (${rmg} ${BLAS_LIBRARIES})
            target_link_libraries (${rmg-tddft} ${BLAS_LIBRARIES})
            target_link_libraries (rmg-on ${BLAS_LIBRARIES})
            target_link_libraries (on-tddft ${BLAS_LIBRARIES})
            target_link_libraries (${rmg-negf} ${BLAS_LIBRARIES})
        endif (BLAS_FOUND)

        find_package(LAPACK REQUIRED)
        if(LAPACK_FOUND)
            target_link_libraries (${rmg} ${LAPACK_LIBRARIES})
            target_link_libraries (${rmg-tddft} ${LAPACK_LIBRARIES})
            target_link_libraries (rmg-on ${LAPACK_LIBRARIES})
            target_link_libraries (on-tddft ${LAPACK_LIBRARIES})
            target_link_libraries (${rmg-negf} ${LAPACK_LIBRARIES})
        endif (LAPACK_FOUND)

        find_package(BLACS)
        if(BLACS_FOUND)
            target_link_libraries (${rmg} ${BLACS_LIBRARIES})
            target_link_libraries (${rmg} ${BLACSCINIT_LIBRARIES})

            target_link_libraries (${rmg-tddft} ${BLACS_LIBRARIES})
            target_link_libraries (${rmg-tddft} ${BLACSCINIT_LIBRARIES})

            target_link_libraries (rmg-on ${BLACS_LIBRARIES})
            target_link_libraries (rmg-on ${BLACSCINIT_LIBRARIES})

            target_link_libraries (on-tddft ${BLACS_LIBRARIES})
            target_link_libraries (on-tddft ${BLACSCINIT_LIBRARIES})
            target_link_libraries (${rmg-negf} ${BLACS_LIBRARIES})     
            target_link_libraries (${rmg-negf} ${BLACSCINIT_LIBRARIES})
        endif(BLACS_FOUND)

        set(USE_SCALAPACK_LIBS 0)
        find_package(SCALAPACK)
        if(SCALAPACK_FOUND)
            unset (USE_SCALAPACK_LIBS CACHE)
            set(USE_SCALAPACK_LIBS 1)
            target_link_libraries (${rmg} ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (${rmg-tddft} ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (rmg-on ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (on-tddft ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (${rmg-negf} ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
        endif (SCALAPACK_FOUND)

    else (NOT CRAY_XK)
        set(USE_SCALAPACK_LIBS 1)
    endif(NOT CRAY_XK)

    # Verify that we can compile and build a scalapack app if the user has set the manual scalapack flag
    if(USE_SCALAPACK_LIBS AND NOT CRAY_XK)
        include("PlatformChecks/CheckScalapack.inc")
        if(NOT RMG_SCALAPACK_INTERNAL)
            unset (USE_SCALAPACK_LIBS CACHE)
            set(USE_SCALAPACK_LIBS 0)
        endif(NOT RMG_SCALAPACK_INTERNAL)
    endif(USE_SCALAPACK_LIBS AND NOT CRAY_XK)

endif(NOT RMG_MATH_LIBS)

if(RMG_PFFT_LIBS)
find_package(PFFT)
    if(NOT PFFT_FOUND)
        set($ENV{MPICC} ${MPI_C_COMPILER})
        execute_process(
        COMMAND tar xvf ${PROJECT_SOURCE_DIR}/lib/pfft-master.tar.gz -C ${PROJECT_SOURCE_DIR}/lib )
        execute_process(
        COMMAND mkdir ${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64)
        execute_process(COMMAND bootstrap.sh 
            WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/lib/pfft-master)
        execute_process(COMMAND
        ${PROJECT_SOURCE_DIR}/lib/pfft-master/configure
        --host=x86_64
        --prefix=${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64
        --libdir=${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64/lib
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64)

        execute_process(COMMAND make  -j 16
            WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/lib/pfft-master)
        execute_process(COMMAND make install
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64)
    endif(NOT PFFT_FOUND)

    find_package(PFFT REQUIRED)
    include_directories(${PFFT_INCLUDES})
    add_definitions(-DUSE_PFFT)
endif(RMG_PFFT_LIBS)


# Look for GPU libraries. Only CUDA supported for now
if(RMG_GPU_ENABLED)
  find_package(CUDA REQUIRED)
  if(${CUDA_VERSION_MAJOR} LESS 6)
      message(FATAL_ERROR "\nCompiling RMG with GPU support requires CUDA version 6.0 or greater. Either set GPU_ENABLED=0 to build a CPU only version or upgrade your cuda installation.") 
  endif(${CUDA_VERSION_MAJOR} LESS 6)

  if(NOT CRAY_XK)
      target_link_libraries(${rmg} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(${rmg-tddft} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(rmg-on ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(on-tddft ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(${rmg-negf} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
  endif(NOT CRAY_XK)
  target_link_libraries(${rmg} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY})
  target_link_libraries(${rmg-tddft} ${CUDA_cublas_LIBRARY})
  target_link_libraries(rmg-on ${CUDA_CUDART_LIBRARY})
  target_link_libraries(rmg-on ${CUDA_cublas_LIBRARY})
  target_link_libraries(on-tddft ${CUDA_CUDART_LIBRARY})
  target_link_libraries(on-tddft ${CUDA_cublas_LIBRARY})
  target_link_libraries(${rmg-negf} ${CUDA_CUDART_LIBRARY})
  target_link_libraries(${rmg-negf} ${CUDA_cublas_LIBRARY})
  include_directories(${CUDA_INCLUDE_DIRS})


  set(CUDA_NVCC_FLAGS "-arch=sm_13" CACHE STRING "nvcc flags" FORCE)
 
#HW capability 1.3
  set (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)
#  set(LIB_TYPE STATIC) #set the lib type
#  CUDA_ADD_LIBRARY(TestLib ${LIB_TYPE}
#  ${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/genvpsi_gpu.cu 
#  ${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/rho_psi_times_psi.cu)
#  target_link_libraries(rmg-negf TestLib)
#  ADD_DEPENDENCIES(rmg-negf TestLib)
 
endif(RMG_GPU_ENABLED)

#include_directories(rmg "/usr/local/magma/include")

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


find_package(LIBSPG)
if(NOT LIBSPG_FOUND)
    execute_process(
    COMMAND tar xvf ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3.tar.gz -C ${PROJECT_SOURCE_DIR}/lib )
    execute_process(
    COMMAND mkdir ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64)
    execute_process(
    COMMAND ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/configure
        --prefix=${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64
        --libdir=${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64)
    execute_process(COMMAND make 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64)
    execute_process(
    	COMMAND cp src/.libs/libsymspg.a . 
    	COMMAND cp src/.libs/libsymspg.la . 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64)
endif(NOT LIBSPG_FOUND)

find_package(LIBSPG REQUIRED)

if(USE_ELEMENTAL_LIBS)
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
    target_link_libraries (rmg-on ON_Common ${ELEMENTAL_LIBRARIES})
    find_package(PMRRR)
    include_directories(${PMRRR_INCLUDES})
    target_link_libraries (rmg-on ${PMRRR_LIBRARIES})
    target_link_libraries (rmg-on ON_Common ${PMRRR_LIBRARIES})
endif(USE_ELEMENTAL_LIBS)


# include files
include_directories(${Boost_INCLUDE_DIRS})
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_SOURCE_DIR}/RmgLib/")
include_directories("${RMG_SOURCE_DIR}/RmgLib/include")
  
include_directories("${PROJECT_SOURCE_DIR}/ON/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/TDDFT/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/NEGF/Headers/")


# Common code subdirectories
add_subdirectory (RmgLib) 
add_subdirectory (Finite_diff) 
add_subdirectory (Force) 
add_subdirectory (Input) 
add_subdirectory (MG) 
add_subdirectory (Misc) 
add_subdirectory (US_PP) 
add_subdirectory (XC_useLIBXC) 
add_subdirectory (XC) 
add_subdirectory (InternalPseudo) 

# Target specific code subdirectories
add_subdirectory (RMG/Common EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin/XC_notuseLIBXC EXCLUDE_FROM_ALL) 
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
link_directories (${RMG_BINARY_DIR}/US_PP)
link_directories (${RMG_BINARY_DIR}/XC)

link_directories (${RMG_BINARY_DIR}/RMG/Common)
link_directories (${RMG_BINARY_DIR}/RMG/Input)
link_directories (${RMG_BINARY_DIR}/RMG/Spin)
link_directories (${RMG_BINARY_DIR}/RMG/Spin/XC)
link_directories (${RMG_BINARY_DIR}/RMG/Subdiag)


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

# Send the configure file to the build directory
configure_file (
  "${PROJECT_SOURCE_DIR}/make_conf.h.in"
  "${PROJECT_BINARY_DIR}/Headers/make_conf.h"
  )

# Build date and time revisions
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    execute_process (COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE foo1) 
    execute_process (COMMAND "cmd" " /C time /T" OUTPUT_VARIABLE foo2) 
    string(CONCAT foo3 ${foo1} ${foo2})
    string(REPLACE "\n" "" BuildDate ${foo3})
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    execute_process (COMMAND date OUTPUT_VARIABLE foo5) 
    string(REPLACE "\n" "" BuildDate ${foo5})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

include(FindSubversion)
if(Subversion_FOUND AND USE_SVN)
    execute_process (
        COMMAND svn info OUTPUT_VARIABLE foo)
    string(REPLACE "Last Changed" ";" foo2 ${foo})
    list(GET foo2 2 foo3)
    string(REPLACE "\n" ";" foo4 ${foo3})
    string(REPLACE " " ";" foo3 ${foo4})
    list(GET foo3 2 CodeRevision)

    list(GET foo2 3 foo3)
    string(REPLACE " " ";" foo4 ${foo3})
    list(GET foo4 2 LastChangeDate)
else(Subversion_FOUND AND USE_SVN)
    set(CodeRevision "release 1.1")
    set(LastChangeDate "Unknown")
endif(Subversion_FOUND AND USE_SVN)


configure_file (
  "${PROJECT_SOURCE_DIR}/svnrev.h.in"
  "${PROJECT_BINARY_DIR}/Headers/svnrev.h"
  )

# These include files include some subproject specific stuff
# in order to keep from cluttering up the main CMakeLists.txt
# file.
include("SubprojectIncludes/CMakeLists.txt.rmg.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-tddft.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-on.inc")
include("SubprojectIncludes/CMakeLists.txt.rmg-negf.inc")
include("SubprojectIncludes/CMakeLists.txt.on-tddft.inc")
if(USE_LINUX_MKL_LIBS)
    target_link_libraries (${rmg-tddft} -Wl,--start-group ${MKL_LIB64}/libmkl_scalapack_lp64.a ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
    target_link_libraries (${rmg} -Wl,--start-group ${MKL_LIB64}/libmkl_scalapack_lp64.a ${MKL_LIB64}/libmkl_intel_lp64.a ${MKL_LIB64}/libmkl_core.a ${MKL_LIB64}/libmkl_gnu_thread.a -Wl,--end-group)
endif(USE_LINUX_MKL_LIBS)

install(TARGETS ${rmg} 
DESTINATION rmg
RUNTIME DESTINATION rmg/bin)

install(FILES RmgLib/libRmgLibShared.so RmgLib/libRmgLib.a
DESTINATION rmg/lib)

# Packaging section next
SET(cpack_input_file "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
SET(CPACK_INPUT_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
SET(CPACK_PROJECT_CONFIG_FILE "${PROJECT_SOURCE_DIR}/CMakeCPackOptions.cmake.in")
set(CPACK_GENERATOR "RPM")

set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(COMPLETE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})
set(CPACK_PACKAGE_NAME "rmg")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Emil Briggs")
set(CPACK_PACKAGE_VENDOR "North Carolina State University")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

include(CPack)

