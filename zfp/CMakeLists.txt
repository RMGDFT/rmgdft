add_library(zfp
src/zfp.c
src/encode4l.c
src/encode4i.c
src/encode4f.c
src/encode4d.c
src/encode3l.c
src/encode3i.c
src/encode3f.c
src/encode3d.c
src/encode2l.c
src/encode2i.c
src/encode2f.c
src/encode2d.c
src/encode1l.c
src/encode1i.c
src/encode1f.c
src/encode1d.c
src/decode4l.c
src/decode4i.c
src/decode4f.c
src/decode4d.c
src/decode3l.c
src/decode3i.c
src/decode3f.c
src/decode3d.c
src/decode2l.c
src/decode2i.c
src/decode2f.c
src/decode2d.c
src/decode1l.c
src/decode1i.c
src/decode1f.c
src/decode1d.c
src/bitstream.c

)


#------------------------------------------------------------------------------#
# Parse version number from zfp.h
#------------------------------------------------------------------------------#
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/include/zfp.h _zfp_h_contents)
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_MAJOR[ \t]+([0-9]+).*"
     "\\1" ZFP_VERSION_MAJOR ${_zfp_h_contents})
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_MINOR[ \t]+([0-9]+).*"
    "\\1" ZFP_VERSION_MINOR ${_zfp_h_contents})
string(REGEX REPLACE ".*#define[ \t]+ZFP_VERSION_PATCH[ \t]+([0-9]+).*"
    "\\1" ZFP_VERSION_PATCH ${_zfp_h_contents})
set(ZFP_VERSION
  "${ZFP_VERSION_MAJOR}.${ZFP_VERSION_MINOR}.${ZFP_VERSION_PATCH}")

#project(ZFP VERSION ${ZFP_VERSION})

#------------------------------------------------------------------------------#
# Top level options
#------------------------------------------------------------------------------#

# Windows (Visual Studio) specific options
if(MSVC)
  # Use this to get a usable export library when building a DLL on Windows
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

  # Silence extraneous Visual Studio specific warnings
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS /wd4146 /wd4305)
endif()

# Suggest C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

if(MSVC OR MINGW)
  set(CMAKE_C_STANDARD 90)
endif()

message(STATUS "Compiling with C standard: ${CMAKE_C_STANDARD}")

# Suggest C++98
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 98)
endif()
message(STATUS "Compiling with C++ standard: ${CMAKE_CXX_STANDARD}")

include(CMakeDependentOption)

# Compile-time options.

set(ZFP_BIT_STREAM_WORD_SIZE 64 CACHE STRING
  "Use smaller bit stream word type for finer rate granularity")
set_property(CACHE ZFP_BIT_STREAM_WORD_SIZE PROPERTY STRINGS "8;16;32;64")

option(ZFP_WITH_CUDA "Enable CUDA parallel compression" OFF)

option(ZFP_WITH_BIT_STREAM_STRIDED
  "Enable strided access for progressive zfp streams" OFF)

option(ZFP_WITH_ALIGNED_ALLOC "Enable aligned memory allocation" OFF)

option(ZFP_WITH_CACHE_TWOWAY "Use two-way skew-associative cache" OFF)

option(ZFP_WITH_CACHE_FAST_HASH
  "Use a faster but more collision prone hash function" OFF)

option(ZFP_WITH_CACHE_PROFILE "Count cache misses" OFF)

# Handle compile-time macros

if((DEFINED ZFP_INT64) AND (DEFINED ZFP_INT64_SUFFIX))
  list(APPEND zfp_public_defs ZFP_INT64=${ZFP_INT64})
  list(APPEND zfp_public_defs ZFP_INT64_SUFFIX=${ZFP_INT64_SUFFIX})
endif()

if((DEFINED ZFP_UINT64) AND (DEFINED ZFP_UINT64_SUFFIX))
  list(APPEND zfp_public_defs ZFP_UINT64=${ZFP_UINT64})
  list(APPEND zfp_public_defs ZFP_UINT64_SUFFIX=${ZFP_UINT64_SUFFIX})
endif()

# This odd cmake pattern here let's the OpenMP feature be either auto-detected,
# explicitly enabled, or explicitly disabled, instead of just on or off.
if(DEFINED ZFP_WITH_OPENMP)
  option(ZFP_WITH_OPENMP "Enable OpenMP parallel compression"
    ${ZFP_WITH_OPENMP})
  if(ZFP_WITH_OPENMP)
    find_package(OpenMP COMPONENTS C REQUIRED)
  endif()
else()
  find_package(OpenMP COMPONENTS C)
  option(ZFP_WITH_OPENMP "Enable OpenMP parallel compression" ${OPENMP_FOUND})
endif()

# Some compilers don't use explicit libraries on the link line for OpenMP but
# instead need to treat the OpenMP C flags as both compile and link flags
# i.e. -fopenmp for compiling and -lgomp for linking, use -fomp for both
# compiling and linking
if(ZFP_WITH_OPENMP AND NOT OpenMP_C_LIBRARIES)
  set(OpenMP_C_LIBRARIES ${OpenMP_C_FLAGS})
endif()

if(ZFP_WITH_CUDA)
  # use CUDA_BIN_DIR hint
  set(ENV{CUDA_BIN_PATH} ${CUDA_BIN_DIR})
  find_package(CUDA)
  if(NOT CUDA_FOUND)
    message(FATAL_ERROR "ZFP_WITH_CUDA is enabled, but a CUDA installation was not found.")
  endif()
  if(${CUDA_VERSION_MAJOR} LESS 7)
        message(FATAL_ERROR "zfp requires at least CUDA 7.0.")
  endif()
endif()

if(NOT (ZFP_BIT_STREAM_WORD_SIZE EQUAL 64))
  list(APPEND zfp_private_defs BIT_STREAM_WORD_TYPE=uint${ZFP_BIT_STREAM_WORD_SIZE})
endif()

if(ZFP_WITH_BIT_STREAM_STRIDED)
  list(APPEND zfp_public_defs BIT_STREAM_STRIDED)
endif()

if(ZFP_WITH_ALIGNED_ALLOC)
  list(APPEND zfp_compressed_array_defs ZFP_WITH_ALIGNED_ALLOC)
endif()

if(ZFP_WITH_CACHE_TWOWAY)
  list(APPEND zfp_compressed_array_defs ZFP_WITH_CACHE_TWOWAY)
endif()

if(ZFP_WITH_CACHE_FAST_HASH)
  list(APPEND zfp_compressed_array_defs ZFP_WITH_CACHE_FAST_HASH)
endif()

if(ZFP_WITH_CACHE_PROFILE)
  list(APPEND zfp_compressed_array_defs ZFP_WITH_CACHE_PROFILE)
endif()

# Link libm only if necessary
include(CheckCSourceCompiles)
check_c_source_compiles("#include<math.h>\nfloat f; int main(){sqrt(f);return 0;}" HAVE_MATH)
if(NOT HAVE_MATH)
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_c_source_compiles("#include<math.h>\nfloat f; int main(){sqrt(f);return 0;}" HAVE_LIBM_MATH)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT HAVE_LIBM_MATH)
    message(FATAL_ERROR "Unable to use C math library functions (with or without -lm)")
  endif()
endif()

#------------------------------------------------------------------------------#
# Add source code
#------------------------------------------------------------------------------#
#include(CTest)
#if(BUILD_TESTING)
#  enable_testing()
#endif()

#set(ZFP_LIBRARY_PREFIX "" CACHE STRING
#  "Prefix to prepend to the output library name")
#mark_as_advanced(ZFP_LIBRARY_PREFIX)

#add_subdirectory(src)

option(BUILD_ALL "Build all subdirectories" OFF)
if(BUILD_ALL)
  set(BUILD_CFP ON CACHE BOOL "Build CFP arrays library" FORCE)
  set(BUILD_ZFORP ON CACHE BOOL "Build Fortran library" FORCE)
  set(BUILD_ZFPY ON CACHE BOOL "Build python bindings for zfp" FORCE)
  set(BUILD_UTILITIES ON CACHE BOOL "Build command line utilities for zfp" FORCE)
  set(BUILD_EXAMPLES ON CACHE BOOL "Build Examples" FORCE)
endif()

option(BUILD_CFP "Build CFP arrays library" OFF)
if(BUILD_CFP)
  add_subdirectory(cfp)
endif()

option(BUILD_ZFORP "Build Fortran library" OFF)
if(BUILD_ZFORP)
  add_subdirectory(fortran)
endif()

option(BUILD_ZFPY "Build python bindings for zfp" OFF)
if(BUILD_ZFPY)
  add_subdirectory(python)
endif()

option(BUILD_UTILITIES "Build command line utilities for zfp" ON)
if(BUILD_UTILITIES)
  add_subdirectory(utils)
endif()

option(BUILD_EXAMPLES "Build Examples" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------#
# Header install
#------------------------------------------------------------------------------#
#install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(DIRECTORY array/   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#if(BUILD_CFP)
#  install(DIRECTORY cfp/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#endif()
#------------------------------------------------------------------------------#
# Build type: one of None, Debug, Release, RelWithDebInfo, MinSizeRel
#------------------------------------------------------------------------------#
#if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE Release)
#endif()

