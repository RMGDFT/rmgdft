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
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required (VERSION 2.6) 
enable_language(Fortran)
project (RMG C CXX Fortran) 
# Load custom cmake modules for FFTW and Magma
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Detect system type
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    set(OperatingSystem "Linux")
    set(RMG_SYSTEM_COMPILER_FLAG "-DLINUX=1")
    set(ZLIB "z")
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

if(DEFINED ENV{CRAY_LIBSCI})
    set (CRAY_XK 1)
endif()


# Compiler flags
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_C_FLAGS "-Wall -Werror-implicit-function-declaration -O3 -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
    set(CMAKE_C_FLAGS "-Wall -O3 ${RMG_SYSTEM_COMPILER_FLAG} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs -O3 -std=c++11 -Wno-write-strings ${RMG_SYSTEM_COMPILER_FLAG} -pthread -D_REENTRANT -DBATCH_NLS=1 -march=native -ffast-math -mfpmath=sse -fopenmp -m64 -I$(MKLROOT)/include")        ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
    set(CMAKE_EXE_LINKER_FLAGS "-ldl -fopenmp")  ## Strip binary
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -march=native -lg2c")
endif()

if (CRAY_XK)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lhugetlbfs")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lhugetlbfs")
endif(CRAY_XK)

set (RMG_GPU_ENABLED 1)
set (RMG_MAGMA_LIBS 1) 

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

set(CMAKE_FIND_ROOT_PATH 
/opt/fftw/default/x86_64     #fftw3 path for bluewater
/usr/local/apps/fftw/intel101_hydra/2.1.5
/usr/local/apps/mvapich2/220b-intel140
/sw/xe/magma/1.5.0/magma-1.5.0
/usr/lib64/mpi/gcc/openmpi/
/sw/xc30/boost/1.55.0/cle5.1_gnu4.8.2
/opt/local/lib/libxc         #libxc from Macports
/home/luw/Downloads/magma-1.4.1)

# Executables defined by the file that contains main()
# Main codes
add_executable(rmg RMG/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-on ON/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-tddft TDDFT/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code
add_executable(rmg-negf NEGF/Common/Main.cpp ${PROJECT_SOURCE_DIR}/Headers/main.h)		# Base code

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


set(USE_PLPLOT_LIBS 0)
find_package(PLPLOT)
if(PLPLOT_FOUND)
    unset (USE_PLPLOT_LIBS CACHE)
    set(USE_PLPLOT_LIBS 1)
    target_link_libraries (rmg ${PLplot_cxx_LIBRARY} ${PLplot_LIBRARY} )
    target_link_libraries (rmg-on ${PLplot_cxx_LIBRARY} ${PLplot_LIBRARY} )
    target_link_libraries (rmg-tddft ${PLplot_cxx_LIBRARY} ${PLplot_LIBRARY} )
    target_link_libraries (rmg-negf ${PLplot_cxx_LIBRARY} ${PLplot_LIBRARY} )
    target_link_libraries (rmg Misc ${PLplot_LIBRARY})
    target_link_libraries (rmg-on Misc ${PLplot_LIBRARY})
    target_link_libraries (rmg-negf Misc ${PLplot_LIBRARY})
    target_link_libraries (rmg-tddft Misc ${PLplot_LIBRARY})
endif(PLPLOT_FOUND)


# static libraries are created for the object files in each subdir. Duplicates
# are needed to handle circular references
target_link_libraries (rmg RmgLib)
target_link_libraries (rmg InternalPseudo)
target_link_libraries (rmg Finite_diff RMG_Common RmgLib ${ZLIB} ${FFTW_LIBRARIES})
target_link_libraries (rmg Misc RMG_Common Finite_diff)
target_link_libraries (rmg RMG_Common Force)
target_link_libraries (rmg Force RMG_Common)
target_link_libraries (rmg Input InternalPseudo)
target_link_libraries (rmg MG RMG_Common RmgLib)
target_link_libraries (rmg Finite_diff Misc)
target_link_libraries (rmg Misc Finite_diff)
target_link_libraries (rmg Input)
target_link_libraries (rmg XC)
target_link_libraries (rmg XC Misc)
target_link_libraries (rmg XC Finite_diff)
target_link_libraries (rmg RMG_Spin XC)
target_link_libraries (rmg XC RMG_Spin)
target_link_libraries (rmg RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg US_PP Input)
target_link_libraries (rmg RMG_Common US_PP)
target_link_libraries (rmg US_PP RMG_Common)
target_link_libraries (rmg Force)
target_link_libraries (rmg Force US_PP)
target_link_libraries (rmg US_PP Force)
target_link_libraries (rmg RMG_Subdiag Finite_diff RMG_Common)
target_link_libraries (rmg RMG_Subdiag Misc)


target_link_libraries (rmg-on RmgLib)
target_link_libraries (rmg-on Finite_diff ON_Common RmgLib)
target_link_libraries (rmg-on Misc ON_Common Finite_diff ${ZLIB} ${FFTW_LIBRARIES})
target_link_libraries (rmg-on ON_Common ON_NEGF_share)
target_link_libraries (rmg-on Misc ON_NEGF_share)
target_link_libraries (rmg-on ON_Common Force)
target_link_libraries (rmg-on Force ON_Common)
target_link_libraries (rmg-on Input ON_Input ON_NEGF_share)
target_link_libraries (rmg-on MG ON_Common RmgLib)
target_link_libraries (rmg-on Finite_diff Misc)
target_link_libraries (rmg-on Misc Finite_diff)
target_link_libraries (rmg-on ON_Input Input)
target_link_libraries (rmg-on XC)
target_link_libraries (rmg-on XC Misc)
target_link_libraries (rmg-on XC Finite_diff)
target_link_libraries (rmg-on RMG_Spin_XC Finite_diff)
target_link_libraries (rmg-on RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg-on RMG_Spin_XC RMG_Spin)
target_link_libraries (rmg-on US_PP InternalPseudo)
target_link_libraries (rmg-on ON_Common)
target_link_libraries (rmg-on Force)
target_link_libraries (rmg-on Force US_PP)
target_link_libraries (rmg-on US_PP Force)
target_link_libraries (rmg-on US_PP Misc)

target_link_libraries (rmg-negf RmgLib)
target_link_libraries (rmg-negf Finite_diff NEGF_Common RmgLib)
target_link_libraries (rmg-negf Misc NEGF_Common Finite_diff ${ZLIB} ${FFTW_LIBRARIES})
target_link_libraries (rmg-negf NEGF_Common ON_NEGF_share)
target_link_libraries (rmg-negf Misc ON_NEGF_share)
target_link_libraries (rmg-negf NEGF_Common Force)
target_link_libraries (rmg-negf Force NEGF_Common)
target_link_libraries (rmg-negf Input NEGF_Input ON_NEGF_share)
target_link_libraries (rmg-negf MG NEGF_Common RmgLib)
target_link_libraries (rmg-negf Finite_diff Misc)
target_link_libraries (rmg-negf Misc Finite_diff)
target_link_libraries (rmg-negf NEGF_Input Input)
target_link_libraries (rmg-negf XC)
target_link_libraries (rmg-negf XC Misc)
target_link_libraries (rmg-negf RMG_Spin_XC)
target_link_libraries (rmg-negf US_PP InternalPseudo)
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
target_link_libraries (rmg-negf NEGF_Pois Finite_diff RmgLib)
target_link_libraries (rmg-negf Misc NEGF_Pois)
target_link_libraries (rmg-negf NEGF_Cudafunc)
target_link_libraries (rmg-negf NEGF_Common NEGF_Cudafunc)
target_link_libraries (rmg-negf NEGF_Cudafunc NEGF_Common)
#target_link_libraries (rmg-negf NEGF_Cudafunc NEGF_Common /home/luw/Downloads/magma-1.4.1/lib/libmagma.a)

target_link_libraries (rmg-tddft RmgLib)
target_link_libraries (rmg-tddft Finite_diff ON_Common RmgLib)
target_link_libraries (rmg-tddft Misc ON_Common Finite_diff ${ZLIB} ${FFTW_LIBRARIES})
target_link_libraries (rmg-tddft ON_Common ON_NEGF_share)
target_link_libraries (rmg-tddft Misc ON_NEGF_share)
target_link_libraries (rmg-tddft ON_Common Force)
target_link_libraries (rmg-tddft Force ON_Common)
target_link_libraries (rmg-tddft Input ON_Input ON_NEGF_share)
target_link_libraries (rmg-tddft MG ON_Common RmgLib)
target_link_libraries (rmg-tddft US_PP InternalPseudo)
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
target_link_libraries (rmg-tddft RMG_Spin_XC Finite_diff)
target_link_libraries (rmg-tddft RMG_Spin RMG_Spin_XC)
target_link_libraries (rmg-tddft RMG_Spin_XC RMG_Spin)


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


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    #message("CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
    set(Boost_USE_STATIC_LIBS        ON)
    set(Boost_USE_MULTITHREADED      ON)
    set(Boost_USE_STATIC_RUNTIME    OFF)
    set(Boost_COMPILER "-vc120")
    set(BOOST_ROOT "C:/local/boost_1_57_0/")
    set(Boost_DEBUG ON)
    set(BOOST_LIBRARYDIR "C:/local/boost_1_57_0/lib64-msvc-12.0/")
    set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0")
    find_package(Boost COMPONENTS thread system iostreams program_options filesystem date_time chrono zlib bzip2 REQUIRED)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    find_package(Boost REQUIRED COMPONENTS thread system iostreams program_options filesystem REQUIRED)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries (rmg ${Boost_LIBRARIES})
target_link_libraries (rmg-on ${Boost_LIBRARIES})
target_link_libraries (rmg-tddft ${Boost_LIBRARIES})
target_link_libraries (rmg-negf ${Boost_LIBRARIES})

# Look for openbabel
set (USE_OPENBABEL_LIBS 0)
find_package(OPENBABEL)
if(OPENBABEL_FOUND)
   target_link_libraries (rmg ${OPENBABEL_LIBRARIES})
   target_link_libraries (rmg-on ${OPENBABEL_LIBRARIES})
   target_link_libraries (rmg-tddft ${OPENBABEL_LIBRARIES})
   target_link_libraries (rmg-negf ${OPENBABEL_LIBRARIES})
   unset (USE_OPENBABEL_LIBS CACHE)
   set(USE_OPENBABEL_LIBS 1)
endif (OPENBABEL_FOUND)

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
        set(MKL_LIBDIR "${MKLROOT}/lib/intel64")
        set(MKL_LIB64 "${MKLROOT}/lib/intel64")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            if(EXISTS "${MKL_LIB64}")
                target_link_libraries (rmg "${MKL_LIB64}/mkl_intel_lp64.lib")
                target_link_libraries (rmg "${MKL_LIB64}/mkl_core.lib")
                target_link_libraries (rmg "${MKL_LIB64}/mkl_intel_thread.lib")
                set(RMG_MATH_LIBS "done")
            else(EXISTS "${MKL_LIB64}")
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif(EXISTS "${MKL_LIB64}")
        else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            if(EXISTS "${MKL_LIB64}")
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_scalapack_lp64.a" "${MKL_LIB64}/libmkl_intel_lp64.a" ${MAGMA_LIBRARIES} )
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_cdft_core.a")
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_intel_lp64.a")
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_core.a")
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_gnu_thread.a" "${MKL_LIB64}/libmkl_intel_lp64.a")
                target_link_libraries (rmg "${MKL_LIB64}/libmkl_blacs_openmpi_lp64.a")
                set(RMG_MATH_LIBS "done")
            else(EXISTS "${MKL_LIB64}")
                message( SEND_ERROR "RMG requires 64 bit MKL." )
            endif(EXISTS "${MKL_LIB64}")
        endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    endif(MKLROOT)
endif(NOT CRAY_XK)


set(USE_SCALAPACK_LIBS 1)
# MKL includes fftw blas, blacs, lapack and scalapack. If it is not present we have to look for these libraries
# individually
if(NOT RMG_MATH_LIBS)

    if(NOT CRAY_XK)

        # FFTW
        find_package(FFTW REQUIRED)
        include_directories(${FFTW_INCLUDES})
        target_link_libraries (rmg ${FFTW_LIBRARIES})
        target_link_libraries (rmg-on ${FFTW_LIBRARIES})
        target_link_libraries (rmg-tddft ${FFTW_LIBRARIES})
        target_link_libraries (rmg-negf ${FFTW_LIBRARIES})
        target_link_libraries (rmg RMG_Common ${FFTW_LIBRARIES})

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
if(RMG_GPU_ENABLED)
  find_package(CUDA REQUIRED)
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


  set(CUDA_NVCC_FLAGS "-arch=sm_13" CACHE STRING "nvcc flags" FORCE)
 
#HW capability 1.3
  set (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)
  set(LIB_TYPE STATIC) #set the lib type
  CUDA_ADD_LIBRARY(TestLib ${LIB_TYPE}
  ${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/genvpsi_gpu.cu 
  ${PROJECT_SOURCE_DIR}/NEGF/Cuda_funcs/rho_psi_times_psi.cu)
  target_link_libraries(rmg-negf TestLib)
  ADD_DEPENDENCIES(rmg-negf TestLib)
 
endif(RMG_GPU_ENABLED)

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
include_directories("${LIBXC_INCLUDES}"  "${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/")
target_link_libraries (rmg ${LIBXC_LIBRARIES})
target_link_libraries (rmg-on ${LIBXC_LIBRARIES})
target_link_libraries (rmg-tddft ${LIBXC_LIBRARIES})
target_link_libraries (rmg-negf ${LIBXC_LIBRARIES})

add_dependencies(rmg libspg)
add_dependencies(rmg-on libspg)
add_dependencies(rmg-tddft libspg)
add_dependencies(rmg-negf libspg)

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
target_link_libraries (rmg ${LIBSPG_LIBRARIES})
target_link_libraries (rmg-on ${LIBSPG_LIBRARIES})
target_link_libraries (rmg-tddft ${LIBSPG_LIBRARIES})
target_link_libraries (rmg-negf ${LIBSPG_LIBRARIES})




# include files
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
if(MAGMA_LIBS)
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
endif(MAGMA_LIBS)


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
if(Subversion_FOUND)
  execute_process (
      COMMAND svn info OUTPUT_VARIABLE foo)
  string(REPLACE "Last Changed" ";" foo2 ${foo})
  list(GET foo2 2 foo3)
  string(REPLACE ":" ";" foo4 ${foo3})
  list(GET foo4 1 CodeRevision)

  list(GET foo2 3 foo3)
  string(REPLACE " " ";" foo4 ${foo3})
  list(GET foo4 2 LastChangeDate)
else(Subversion_FOUND)
  set(CodeRevision "\"release 1.1\"")
endif(Subversion_FOUND)


configure_file (
  "${PROJECT_SOURCE_DIR}/svnrev.h.in"
  "${PROJECT_BINARY_DIR}/Headers/svnrev.h"
  )

install(TARGETS rmg 
DESTINATION rmg
RUNTIME DESTINATION rmg/bin)

install(FILES RmgLib/libRmgLibShared.so RmgLib/libRmgLib.a
DESTINATION rmg/lib)

# Packaging section next
set(CPACK_GENERATOR "DEB")

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(COMPLETE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})
set(CPACK_PACKAGE_NAME "rmg")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Emil Briggs")
set(CPACK_PACKAGE_VENDOR "North Carolina State University")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

if(CPACK_GENERATOR STREQUAL "RPM")
# RPM specific tags
    set(CPACK_RPM_PACKAGE_RELOCATABLE "ON")
    set(CPACK_RPM_PACKAGE_SUMMARY "Ab-initio dft based electronic structure code")
    set(CPACK_RPM_PACKAGE_DESCRIPTION "RMG is an ab-initio dft based electronic structure code that uses real space
 grids to represent quantities of interest such as the electronic orbitals and
charge density.")
    set(CPACK_RPM_PACKAGE_LICENSE "GPLv2")
    set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
    set(CPACK_RPM_PACKAGE_RELEASE "OpenSuse_13.1")
#set(CPACK_RPM_PACKAGE_REQUIRES "liblapack3 >= 3.4.0,libopenblas0 >= 0.2.8")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${COMPLETE_VERSION}-${CPACK_PACKAGE_RELEASE}-${CPACK_RPM_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

endif(CPACK_GENERATOR STREQUAL "RPM")

if(CPACK_GENERATOR STREQUAL "DEB")

    set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "RMG is an ab-initio dft based electronic structure code that uses real space grids to represent quantities of interest such as the electronic orbitals and charge density.")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libplplot12, libfftw3-double3, libmpich10, libboost-thread1.55.0, libboost-system1.55.0, libboost-iostreams1.55.0, libboost-program-options1.55.0, libboost-filesystem1.55.0, libopenbabel4, libopenblas-base, libgfortran3, libstdc++6, libgomp1, libgcc1, libc6, plplot12-driver-cairo")

endif(CPACK_GENERATOR STREQUAL "DEB")

include(CPack)

