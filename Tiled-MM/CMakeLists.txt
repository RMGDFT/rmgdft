cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

# generate compile_commands.json with compile commands for each target
set(CMAKE_EXPORT_COMPILE_COMMANDS "YES")

set(MASTER_PROJECT OFF)
if(NOT DEFINED PROJECT_NAME)
    set(MASTER_PROJECT ON)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
project(Tiled-MM VERSION 2.3 LANGUAGES CXX)

option(TILEDMM_WITH_TESTS "Enable tests." ${MASTER_PROJECT})
option(TILEDMM_WITH_EXAMPLES "Enable examples." ${MASTER_PROJECT})
set(TILEDMM_GPU_BACKEND "ROCM" CACHE STRING "GPU backend")
set_property(CACHE TILEDMM_GPU_BACKEND PROPERTY STRINGS "CUDA" "ROCM")

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

if(MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release." FORCE)
endif()

message(STATUS "Selected TILEDMM_GPU_BACKEND: ${TILEDMM_GPU_BACKEND}")

set (TILEDMM_DEFINITIONS)
set (TILEDMM_EXTERNAL_LIBRARIES)
# Dependencies
#
if(TILEDMM_GPU_BACKEND STREQUAL "CUDA")
    find_package(CUDAToolkit REQUIRED)
    # this seems to be a cmake bug that was fixed in 3.12.2
    if(CUDA_VERSION VERSION_GREATER 9.1 AND CMAKE_VERSION VERSION_LESS 3.12.2)
         list(REMOVE_ITEM CUDA_CUBLAS_LIBRARIES "CUDA_cublas_device_LIBRARY-NOTFOUND")
    endif()
    if (NOT TARGET tiledmm::cuda)
        add_library(tiledmm::cuda INTERFACE IMPORTED)
        target_link_libraries(tiledmm::cuda INTERFACE CUDA::cublas CUDA::cudart)
    endif()
    list(APPEND TILEDMM_DEFINITIONS TILED_MM_CUDA)

elseif(TILEDMM_GPU_BACKEND STREQUAL "ROCM")
    find_package(rocblas REQUIRED)
    if (NOT TARGET tiledmm::rocm)
      add_library(tiledmm::rocm INTERFACE IMPORTED)
      target_link_libraries(tiledmm::rocm INTERFACE roc::rocblas)
    endif()

    list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCM)
    # TODO: it should be possible to remove this (and adapt the cpp code
    # accordingly). Recent rocblas supports all types
    list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCBLAS_HAS_SGEMM
      TILED_MM_ROCBLAS_HAS_DGEMM TILED_MM_ROCBLAS_HAS_ZGEMM)
else()
    message(FATAL_ERROR "TILED_MM: GPU Backend must be CUDA or ROCM!")
endif()

# Tiled-MM
#
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

add_subdirectory(src/Tiled-MM)

install(EXPORT Tiled-MM_targets
  FILE Tiled-MMTargets.cmake
  NAMESPACE Tiled-MM::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tiled-MM")

install(DIRECTORY "${Tiled-MM_SOURCE_DIR}/src/Tiled-MM"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.hpp")

write_basic_package_version_file("${Tiled-MM_BINARY_DIR}/Tiled-MMConfigVersion.cmake"
  VERSION ${Tiled-MM_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_file("${Tiled-MM_SOURCE_DIR}/cmake/Tiled-MMConfig.cmake.in"
  "${Tiled-MM_BINARY_DIR}/Tiled-MMConfig.cmake"
  @ONLY)

install(FILES "${Tiled-MM_BINARY_DIR}/Tiled-MMConfig.cmake"
  "${Tiled-MM_BINARY_DIR}/Tiled-MMConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tiled-MM")

if(TILEDMM_WITH_TESTS OR TILEDMM_WITH_EXAMPLES)
  find_package(cxxopts REQUIRED)
endif()

if(TILEDMM_WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(TILEDMM_WITH_EXAMPLES)
  add_subdirectory(examples)
endif()
