#
#    RMG master cmake file
#
#    Subprojects:
#        rmg        Base code
#        rmg-on     ON code
#        NEGF
#        rmg-tddft:  TDDFT based on rmg-on
#
# Environment variables that may also control configuration
# MKLROOT = root directory of Intel MKL libraries if installed.
#
#
#
cmake_minimum_required (VERSION 2.6) 
enable_language(Fortran)
project (RMG C CXX Fortran) 

# Detect system type
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    SET(OperatingSystem "Linux")
    SET(RMG_SYSTEM_COMPILER_FLAG "-DLINUX=1")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(DEFINED ENV{CRAY_LIBSCI})
    set (CRAY_XK 1)
#set (RMG_GPU_ENABLED 1)
endif()

set(CMAKE_FIND_ROOT_PATH 
/opt/fftw/default/x86_64     #fftw3 path for bluewater
/usr/local/apps/fftw/intel101_hydra/2.1.5
/usr/local/apps/mvapich2/220b-intel140
/sw/xe/magma/1.3.0/cnl4.2_gnu4.7.2
/usr/lib64/mpi/gcc/openmpi/
/sw/xc30/boost/1.55.0/cle5.1_gnu4.8.2
/home/luw/Downloads/magma-1.4.1
${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64)

# Executables defined by the file that contains main()
# Main codes
add_executable(rmg RMG/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-on ON/Common/md.c ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-tddft TDDFT/Common/main.c ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-negf NEGF/Common/md.c ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code

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

# Load custom cmake modules for FFTW and Magma
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)




# static libraries are created for the object files in each subdir. Duplicates
# are needed to handle circular references
target_link_libraries (rmg RmgLib)
target_link_libraries (rmg Finite_diff RMG_Common RmgLib)
target_link_libraries (rmg Misc RMG_Common Finite_diff)
target_link_libraries (rmg RMG_Common Force)
target_link_libraries (rmg Force RMG_Common)
target_link_libraries (rmg Input RMG_Input)
target_link_libraries (rmg MG RMG_Common RmgLib)
target_link_libraries (rmg Misc)
target_link_libraries (rmg Finite_diff Misc)
target_link_libraries (rmg Misc Finite_diff)
target_link_libraries (rmg RMG_Input Input)
target_link_libraries (rmg XC)
target_link_libraries (rmg XC Misc)
target_link_libraries (rmg XC Finite_diff)
target_link_libraries (rmg RMG_Spin XC)
target_link_libraries (rmg XC RMG_Spin)
target_link_libraries (rmg RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg US_PP)
target_link_libraries (rmg RMG_Common ${FFTW_LIBRARIES})
target_link_libraries (rmg Force)
target_link_libraries (rmg Force US_PP)
target_link_libraries (rmg US_PP Force)
target_link_libraries (rmg RMG_Subdiag Finite_diff)
target_link_libraries (rmg RMG_Subdiag Misc)


target_link_libraries (rmg-on RmgLib)
target_link_libraries (rmg-on Finite_diff ON_Common RmgLib)
target_link_libraries (rmg-on Misc ON_Common Finite_diff)
target_link_libraries (rmg-on ON_Common ON_NEGF_share)
target_link_libraries (rmg-on Misc ON_NEGF_share)
target_link_libraries (rmg-on ON_Common Force)
target_link_libraries (rmg-on Force ON_Common)
target_link_libraries (rmg-on Input ON_Input ON_NEGF_share)
target_link_libraries (rmg-on MG ON_Common RmgLib)
target_link_libraries (rmg-on Misc)
target_link_libraries (rmg-on Finite_diff Misc)
target_link_libraries (rmg-on Misc Finite_diff)
target_link_libraries (rmg-on ON_Input Input)
target_link_libraries (rmg-on XC)
target_link_libraries (rmg-on XC Misc)
target_link_libraries (rmg-on XC Finite_diff)
target_link_libraries (rmg-on RMG_Spin_XC Finite_diff)
target_link_libraries (rmg-on RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg-on RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg-on US_PP)
target_link_libraries (rmg-on ON_Common)
target_link_libraries (rmg-on Force)
target_link_libraries (rmg-on Force US_PP)
target_link_libraries (rmg-on US_PP Force)
target_link_libraries (rmg-on US_PP Misc)

target_link_libraries (rmg-negf RmgLib)
target_link_libraries (rmg-negf Finite_diff NEGF_Common RmgLib)
target_link_libraries (rmg-negf Misc NEGF_Common Finite_diff)
target_link_libraries (rmg-negf NEGF_Common ON_NEGF_share)
target_link_libraries (rmg-negf Misc ON_NEGF_share)
target_link_libraries (rmg-negf NEGF_Common Force)
target_link_libraries (rmg-negf Force NEGF_Common)
target_link_libraries (rmg-negf Input NEGF_Input ON_NEGF_share)
target_link_libraries (rmg-negf MG NEGF_Common RmgLib)
target_link_libraries (rmg-negf Misc)
target_link_libraries (rmg-negf Finite_diff Misc)
target_link_libraries (rmg-negf Misc Finite_diff)
target_link_libraries (rmg-negf NEGF_Input Input)
target_link_libraries (rmg-negf XC)
target_link_libraries (rmg-negf XC Misc)
target_link_libraries (rmg-negf RMG_Spin_XC)
target_link_libraries (rmg-negf US_PP)
target_link_libraries (rmg-negf NEGF_Common Misc)
target_link_libraries (rmg-negf ON_NEGF_share NEGF_Common)
target_link_libraries (rmg-negf NEGF_Common)
target_link_libraries (rmg-negf ON_NEGF_share)
target_link_libraries (rmg-negf Force)
target_link_libraries (rmg-negf Force US_PP)
target_link_libraries (rmg-negf US_PP Force)
target_link_libraries (rmg-negf US_PP Misc)
target_link_libraries (rmg-negf NEGF_Lcurrent)
target_link_libraries (rmg-negf NEGF_Pois Misc)
target_link_libraries (rmg-negf NEGF_Pois MG RmgLib)
target_link_libraries (rmg-negf NEGF_Pois Finite_diff)
target_link_libraries (rmg-negf Misc NEGF_Pois)
target_link_libraries (rmg-negf NEGF_Cudafunc)
target_link_libraries (rmg-negf NEGF_Common NEGF_Cudafunc)
target_link_libraries (rmg-negf NEGF_Cudafunc NEGF_Common)
#target_link_libraries (rmg-negf NEGF_Cudafunc NEGF_Common /home/luw/Downloads/magma-1.4.1/lib/libmagma.a)

target_link_libraries (rmg-tddft RmgLib)
target_link_libraries (rmg-tddft Finite_diff ON_Common RmgLib)
target_link_libraries (rmg-tddft Misc ON_Common Finite_diff)
target_link_libraries (rmg-tddft ON_Common ON_NEGF_share)
target_link_libraries (rmg-tddft Misc ON_NEGF_share)
target_link_libraries (rmg-tddft ON_Common Force)
target_link_libraries (rmg-tddft Force ON_Common)
target_link_libraries (rmg-tddft Input ON_Input ON_NEGF_share)
target_link_libraries (rmg-tddft MG ON_Common RmgLib)
target_link_libraries (rmg-tddft Misc)
target_link_libraries (rmg-tddft Finite_diff Misc)
target_link_libraries (rmg-tddft Misc Finite_diff)
target_link_libraries (rmg-tddft ON_Input Input)
target_link_libraries (rmg-tddft XC)
target_link_libraries (rmg-tddft XC Misc)
target_link_libraries (rmg-tddft RMG_Spin_XC)
target_link_libraries (rmg-tddft RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg-tddft RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg-tddft US_PP)
target_link_libraries (rmg-tddft ON_Common)
target_link_libraries (rmg-tddft RMG_Spin_XC)
target_link_libraries (rmg-tddft tddft_Common XC RMG_Spin_XC)
target_link_libraries (rmg-tddft RMG_Spin_XC XC)
target_link_libraries (rmg-tddft XC RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg-tddft RMG_Spin_XC RMG_Spin XC)
target_link_libraries (rmg-tddft Force)
target_link_libraries (rmg-tddft Force US_PP)
target_link_libraries (rmg-tddft US_PP Force)
target_link_libraries (rmg-tddft US_PP Misc)
target_link_libraries (rmg-tddft tddft_Common ON_Common Eldyn)
target_link_libraries (rmg-tddft tddft_Common ON_NEGF_share)
target_link_libraries (rmg-tddft tddft_Common Misc)
target_link_libraries (rmg-tddft ON_NEGF_share XC)
target_link_libraries (rmg-tddft Eldyn)
target_link_libraries (rmg-tddft Misc tddft_Common Finite_diff)
target_link_libraries (rmg-tddft tddft_Common ON_NEGF_share)
target_link_libraries (rmg-tddft MG tddft_Common RmgLib)
target_link_libraries (rmg-tddft ON_Common US_PP)
target_link_libraries (rmg-tddft US_PP Misc)
target_link_libraries (rmg-tddft ON_Common XC)


# Search for an MPI implementation
#set(MPI_C_COMPILER "/opt/openmpi/bin/mpicc")
#find_package(MPI REQUIRED)
find_package(MPI)
  include_directories(${MPI_INCLUDE_PATH})
if(MPI_CXX_FOUND)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  target_link_libraries (rmg ${MPI_CXX_LIBRARIES})
  target_link_libraries (rmg-on ${MPI_CXX_LIBRARIES})
  target_link_libraries (rmg-tddft ${MPI_CXX_LIBRARIES})
  target_link_libraries (rmg-negf ${MPI_CXX_LIBRARIES})
else (MPI_CXX_FOUND)
    include(CMakeForceCompiler)
    if(CRAY_XK)
	CMAKE_FORCE_CXX_COMPILER(CC ${CMAKE_CXX_COMPILER_ID})
    endif(CRAY_XK)	
endif (MPI_CXX_FOUND)
if(MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  target_link_libraries (rmg ${MPI_C_LIBRARIES})
  target_link_libraries (rmg-on ${MPI_C_LIBRARIES})
  target_link_libraries (rmg-tddft ${MPI_C_LIBRARIES})
  target_link_libraries (rmg-negf ${MPI_C_LIBRARIES})
else(MPI_C_FOUND)
  CMAKE_FORCE_C_COMPILER(cc ${CMAKE_C_COMPILER_ID})
endif (MPI_C_FOUND)

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
target_link_libraries (rmg ${FFTW_LIBRARIES})
target_link_libraries (rmg-on ${FFTW_LIBRARIES})
target_link_libraries (rmg-tddft ${FFTW_LIBRARIES})
target_link_libraries (rmg-negf ${FFTW_LIBRARIES})

find_package(Boost REQUIRED COMPONENTS thread system)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries (rmg ${Boost_LIBRARIES})
target_link_libraries (rmg-on ${Boost_LIBRARIES})
target_link_libraries (rmg-tddft ${Boost_LIBRARIES})
target_link_libraries (rmg-negf ${Boost_LIBRARIES})


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

    if(NOT CRAY_XK)

        # First look for openblas
        SET(OPEN_BLAS_SEARCH_PATHS  /lib/ /lib64/  /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /opt/OpenBLAS/lib $ENV{OPENBLAS_HOME}/lib )
        FIND_LIBRARY(OPENBLAS NAMES openblas PATHS ${OPEN_BLAS_SEARCH_PATHS})
        IF (OPENBLAS)
            set(BLAS_FOUND TRUE)
            set(BLAS_LIBRARIES ${OPENBLAS})
            set(BLAS_blas_LIBRARY ${OPENBLAS})
        ELSE (OPENBLAS)
            find_package(BLAS REQUIRED)
        ENDIF (OPENBLAS)

        if(BLAS_FOUND)
            target_link_libraries (rmg ${BLAS_LIBRARIES})
            target_link_libraries (rmg-on ${BLAS_LIBRARIES})
            target_link_libraries (rmg-tddft ${BLAS_LIBRARIES})
            target_link_libraries (rmg-negf ${BLAS_LIBRARIES})
        endif (BLAS_FOUND)

        find_package(LAPACK REQUIRED)
        if(LAPACK_FOUND)
            target_link_libraries (rmg ${LAPACK_LIBRARIES})
            target_link_libraries (rmg-on ${LAPACK_LIBRARIES})
            target_link_libraries (rmg-tddft ${LAPACK_LIBRARIES})
            target_link_libraries (rmg-negf ${LAPACK_LIBRARIES})
        endif (LAPACK_FOUND)

        find_package(BLACS)
        if(BLACS_FOUND)
            target_link_libraries (rmg ${BLACS_LIBRARIES})
            target_link_libraries (rmg ${BLACSCINIT_LIBRARIES})

            target_link_libraries (rmg-on ${BLACS_LIBRARIES})
            target_link_libraries (rmg-on ${BLACSCINIT_LIBRARIES})

            target_link_libraries (rmg-tddft ${BLACS_LIBRARIES})
            target_link_libraries (rmg-tddft ${BLACSCINIT_LIBRARIES})
            target_link_libraries (rmg-negf ${BLACS_LIBRARIES})     
            target_link_libraries (rmg-negf ${BLACSCINIT_LIBRARIES})
        endif(BLACS_FOUND)

        set(USE_SCALAPACK_LIBS 0)
        find_package(SCALAPACK)
        if(SCALAPACK_FOUND)
            unset (USE_SCALAPACK_LIBS CACHE)
            set(USE_SCALAPACK_LIBS 1)
            target_link_libraries (rmg ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (rmg-on ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (rmg-tddft ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
            target_link_libraries (rmg-negf ${SCALAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
        endif (SCALAPACK_FOUND)

    else (NOT CRAY_XK)
        set(USE_SCALAPACK_LIBS 1)
    endif(NOT CRAY_XK)

endif(NOT RMG_MATH_LIBS)


# Look for GPU libraries. Only CUDA supported for now
find_package(CUDA)
if(CUDA_FOUND)
  if(NOT CRAY_XK)
      target_link_libraries(rmg ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(rmg-on ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(rmg-tddft ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
      target_link_libraries(rmg-negf ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
  endif(NOT CRAY_XK)
  target_link_libraries(rmg ${CUDA_cublas_LIBRARY})
  target_link_libraries(rmg-on ${CUDA_CUDART_LIBRARY})
  target_link_libraries(rmg-on ${CUDA_cublas_LIBRARY})
  target_link_libraries(rmg-tddft ${CUDA_CUDART_LIBRARY})
  target_link_libraries(rmg-tddft ${CUDA_cublas_LIBRARY})
  target_link_libraries(rmg-negf ${CUDA_CUDART_LIBRARY})
  target_link_libraries(rmg-negf ${CUDA_cublas_LIBRARY})
  include_directories(${CUDA_INCLUDE_DIRS})


SET(CUDA_NVCC_FLAGS "-arch=sm_13" CACHE STRING "nvcc flags" FORCE)
 
#HW capability 1.3
SET (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)
SET(LIB_TYPE STATIC) #set the lib type
CUDA_ADD_LIBRARY(TestLib ${LIB_TYPE}
${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/genvpsi_gpu.cu 
${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/rho_psi_times_psi.cu)
  target_link_libraries(rmg-negf TestLib)
ADD_DEPENDENCIES(rmg-negf TestLib)
 
endif(CUDA_FOUND)

#include_directories(rmg "/usr/local/magma/include")
add_dependencies(rmg libxc)
add_dependencies(rmg-on libxc)
add_dependencies(rmg-tddft libxc)
add_dependencies(rmg-negf libxc)

find_package(LIBXC)
if(NOT LIBXC_FOUND)
    execute_process(
    COMMAND tar xvf ${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3.tar -C ${PROJECT_SOURCE_DIR}/lib )
    execute_process(
    COMMAND mkdir ${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64)
    execute_process(
    COMMAND ${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/configure
        --prefix=${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64
        --libdir=${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64/lib
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64)
    execute_process(COMMAND make install
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64)
endif(NOT LIBXC_FOUND)


find_package(LIBXC REQUIRED)
include_directories("${LIBXC_INCLUDES}")
target_link_libraries (rmg ${LIBXC_LIBRARIES})
target_link_libraries (rmg-on ${LIBXC_LIBRARIES})
target_link_libraries (rmg-tddft ${LIBXC_LIBRARIES})
target_link_libraries (rmg-negf ${LIBXC_LIBRARIES})



# Compiler flags
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_C_FLAGS "-Wall -Werror-implicit-function-declaration -O3 -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -mavx -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
    set(CMAKE_C_FLAGS "-Wall -O3 ${RMG_SYSTEM_COMPILER_FLAG} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -mavx -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs -O3 -std=c++11 ${RMG_SYSTEM_COMPILER_FLAG} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -mavx -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
    set(CMAKE_EXE_LINKER_FLAGS "-ldl -fopenmp")  ## Strip binary
endif()
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c")

if (CRAY_XK)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lhugetlbfs")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lhugetlbfs")
endif(CRAY_XK)

#message(${CMAKE_C_FLAGS})
# include files
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_SOURCE_DIR}/RmgLib/")
include_directories("${RMG_SOURCE_DIR}/RmgLib/include")
  
include_directories("${PROJECT_SOURCE_DIR}/ON/Headers/")
include_directories("${PROJECT_SOURCE_DIR}/NEGF/Headers/")
#include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
#include_directories("${RMG_BINARY_DIR}/RMG/Headers/")
#include_directories("/opt/openmpi/include")


# Common code subdirectories
add_subdirectory (RmgLib) 
add_subdirectory (Finite_diff) 
add_subdirectory (Force) 
add_subdirectory (Input) 
add_subdirectory (MG) 
add_subdirectory (Misc) 
add_subdirectory (US_PP) 
add_subdirectory (XC_useLIBXC) 

# Target specific code subdirectories
add_subdirectory (RMG/Common EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Input EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Spin/XC_notuseLIBXC EXCLUDE_FROM_ALL) 
add_subdirectory (RMG/Subdiag EXCLUDE_FROM_ALL)

add_subdirectory (ON/Common EXCLUDE_FROM_ALL) 
add_subdirectory (ON/Input EXCLUDE_FROM_ALL) 
add_subdirectory (ON/ON-NEGF-share EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/Common EXCLUDE_FROM_ALL) 
add_subdirectory (TDDFT/ELDYN EXCLUDE_FROM_ALL) 

add_subdirectory (NEGF/Common EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Input EXCLUDE_FROM_ALL) 
add_subdirectory (NEGF/Cuda_funcs EXCLUDE_FROM_ALL) 
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

#target_link_libraries (rmg "/usr/local/magma/lib/libmagma.a" "${MKL_LIB64}/libmkl_scalapack_lp64.a")
#target_link_libraries (rmg "/usr/local/magma/lib/libmagma.a" "${MKL_LIB64}/libmkl_intel_lp64.a")
#target_link_libraries (rmg "/usr/local/magma/lib/libmagma.a" "${MKL_LIB64}/libmkl_intel_lp64.a" "/usr/local/magma/lib/libmagma.a")
#target_link_libraries (rmg "/usr/local/magma/lib/libmagma.a" "${MKL_LIB64}/libmkl_gnu_thread.a" "/usr/local/magma/lib/libmagma.a")

# Look for the magma libraries
find_package(MAGMA)
include_directories(${MAGMA_INCLUDES})
if(MAGMA_FOUND)
    if(MAGMABLAS_FOUND) 
        target_link_libraries (rmg ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (rmg-on ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (rmg-tddft ${MAGMA_LIBRARIES}  ${MAGMABLAS_LIBRARIES} )
        target_link_libraries (rmg-negf ${MAGMA_LIBRARIES} ${MAGMABLAS_LIBRARIES} )
    else(MAGMABLAS_FOUND) 
        target_link_libraries (rmg ${MAGMA_LIBRARIES} )
        target_link_libraries (rmg-on ${MAGMA_LIBRARIES} )
        target_link_libraries (rmg-tddft ${MAGMA_LIBRARIES} )
        target_link_libraries (rmg-negf ${MAGMA_LIBRARIES} )
    endif(MAGMABLAS_FOUND) 
endif(MAGMA_FOUND)


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

set (FAST_ORTHO 0)
if(EXPERIMENTAL_FEATURES)
    option (RMG_FAST_ORTHO "Experimental fast ortho. Disabled by default." OFF)
    if (RMG_FAST_ORTHO)
        set (FAST_ORTHO 1)
    else (RMG_FAST_ORTHO)
        set (FAST_ORTHO 0)
    endif(RMG_FAST_ORTHO)
endif(EXPERIMENTAL_FEATURES)

option (RMG_GPU_ENABLED "Use GPU acceleration." ON)
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
    set (MAGMA_LIBS 0)
    if (RMG_GPU_ENABLED)
        set (MAGMA_LIBS 1)
    endif(RMG_GPU_ENABLED)
else (RMG_MAGMA_LIBS)
    set (MAGMA_LIBS 0)
endif(RMG_MAGMA_LIBS)

# Send the configure file to the build directory
configure_file (
  "${PROJECT_SOURCE_DIR}/make_conf.h.in"
  "${PROJECT_BINARY_DIR}/Headers/make_conf.h"
  )

execute_process (
    COMMAND svn info OUTPUT_VARIABLE foo)
string(REPLACE "Last Changed" ";" foo2 ${foo})
list(GET foo2 2 foo3)
string(REPLACE ":" ";" foo4 ${foo3})
list(GET foo4 1 CodeRevision)

list(GET foo2 3 foo3)
string(REPLACE " " ";" foo4 ${foo3})
list(GET foo4 2 LastChangeDate)
execute_process (
    COMMAND date OUTPUT_VARIABLE foo5)
string(REPLACE "\n" "" BuildDate ${foo5})

configure_file (
  "${PROJECT_SOURCE_DIR}/svnrev.h.in"
  "${PROJECT_BINARY_DIR}/Headers/svnrev.h"
  )

