cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

set(MASTER_PROJECT OFF)
if(NOT DEFINED PROJECT_NAME)
    set(MASTER_PROJECT ON)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
project(Tiled-MM VERSION 2.0 LANGUAGES CXX)

option(TILEDMM_WITH_TESTS "Enable tests." ${MASTER_PROJECT})
option(TILEDMM_WITH_EXAMPLES "Enable examples." ${MASTER_PROJECT})
option(TILEDMM_WITH_INSTALL "Enable installation." ${MASTER_PROJECT})
if (NOT TILEDMM_GPU_BACKEND)
    set(TILEDMM_GPU_BACKEND "ROCM" CACHE STRING "GPU backend")
endif()
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
    list(APPEND TILEDMM_EXTERNAL_LIBRARIES tiledmm::cuda)
    list(APPEND TILEDMM_DEFINITIONS TILED_MM_CUDA)

elseif(TILEDMM_GPU_BACKEND STREQUAL "ROCM")
    find_package(HIPLIBS REQUIRED)
    find_package(ROCBLAS REQUIRED)
    list(APPEND TILEDMM_EXTERNAL_LIBRARIES ROCBLAS::rocblas HIPLIBS::hiplibs)
    list(APPEND TILEDMM_DEFINITIONS -D__HIP_PLATFORM_HCC__) # since we are not using the hcc compiler, we have to set up this explicitly
    list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCM)
    if(ROCBLAS_HAS_SGEMM)
        list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCBLAS_HAS_SGEMM)
    endif()
    if(ROCBLAS_HAS_DGEMM)
        list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCBLAS_HAS_DGEMM)
    endif()
    if(ROCBLAS_HAS_CGEMM)
        list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCBLAS_HAS_CGEMM)
    endif()
    if(ROCBLAS_HAS_ZGEMM)
        list(APPEND TILEDMM_DEFINITIONS TILED_MM_ROCBLAS_HAS_ZGEMM)
    endif()
else()
    message(FATAL_ERROR "TILED_MM: GPU Backend must be CUDA or ROCM!")
endif()

# Tiled-MM
#
add_subdirectory(src/Tiled-MM)

if(TILEDMM_WITH_INSTALL)
    include(CMakePackageConfigHelpers)
    include(GNUInstallDirs)

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
                  "${Tiled-MM_SOURCE_DIR}/cmake/FindROCBLAS.cmake"
                  "${Tiled-MM_SOURCE_DIR}/cmake/FindHIPLIBS.cmake"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tiled-MM")
endif()

if(TILEDMM_WITH_TESTS OR TILEDMM_WITH_EXAMPLES)
    add_subdirectory(libs/cxxopts)
endif()

if(TILEDMM_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(TILEDMM_WITH_EXAMPLES)
    add_subdirectory(examples)
endif()
