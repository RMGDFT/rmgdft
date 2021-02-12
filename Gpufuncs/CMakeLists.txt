if(RMG_CUDA_ENABLED)
# This is necessary in order to work around a cmake bug
# in version prior to 3.10 where quotes around string
# literals are escaped in a way that causes nvcc to fail.
remove_definitions(-DBUILD_DATE="${BuildDate}")

# c++11 standard is required for the main code but causes nvcc to fail
# so don't propagate the flags.
set(CUDA_PROPAGATE_HOST_FLAGS OFF )

# Borrowed from the magma libraries cmake files
set(GPU_TARGET "sm_60")
set( NV_SM    "" )
set( NV_COMP  "" )
if ( ${GPU_TARGET} MATCHES sm_60 )
	if ( NOT MIN_ARCH )
		set( MIN_ARCH 600 )
	endif()
	set( NV_SM "${NV_SM} -gencode arch=compute_60,code=sm_60" )
	set( NV_COMP        "-gencode arch=compute_60,code=compute_60" )
	message( STATUS "    compile for CUDA arch 6.0 (Pascal)" )
endif()

if ( ${GPU_TARGET} MATCHES sm_61 )
	if ( NOT MIN_ARCH )
		set( MIN_ARCH 610 )
	endif()
	set( NV_SM "${NV_SM} -gencode arch=compute_61,code=sm_61" )
	set( NV_COMP        "-gencode arch=compute_61,code=compute_61" )
	message( STATUS "    compile for CUDA arch 6.1 (Pascal)" )
endif()

if ( ${GPU_TARGET} MATCHES sm_62 )
	if ( NOT MIN_ARCH )
		set( MIN_ARCH 620 )
	endif()
	set( NV_SM "${NV_SM} -gencode arch=compute_62,code=sm_62" )
	set( NV_COMP        "-gencode arch=compute_62,code=compute_62" )
	message( STATUS "    compile for CUDA arch 6.2 (Pascal)" )
endif()

if ( ${GPU_TARGET} MATCHES sm_70 )
	if ( NOT MIN_ARCH )
		set( MIN_ARCH 700 )
	endif()
	set( NV_SM "${NV_SM} -gencode arch=compute_70,code=sm_70" )
	set( NV_COMP        "-gencode arch=compute_70,code=compute_70" )
	message( STATUS "    compile for CUDA arch 7.0 (Volta)" )
endif()

if ( ${GPU_TARGET} MATCHES sm_71 )
	if ( NOT MIN_ARCH )
		set( MIN_ARCH 710 )
	endif()
	set( NV_SM "${NV_SM} -gencode arch=compute_71,code=sm_71" )
	set( NV_COMP        "-gencode arch=compute_71,code=compute_71" )
	message( STATUS "    compile for CUDA arch 7.1 (Volta)" )
endif()

if ( NOT MIN_ARCH )
        message( FATAL_ERROR "GPU_TARGET must contain one or more of Pascal, Volta, or valid sm_[0-9][0-9]" )
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CUDA_NVCC_FLAGS_RELEASE "${CUDA_NVCC_FLAGS_RELEASE} ${NV_SM} ${NV_COMP} -DCUDA_API_PER_THREAD_DEFAULT_STREAM=1")
elseif (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} ${NV_SM} ${NV_COMP} -DCUDA_API_PER_THREAD_DEFAULT_STREAM=1")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

add_definitions( "-DMIN_CUDA_ARCH=${MIN_ARCH}" )
set(LIB_TYPE STATIC)
  CUDA_ADD_LIBRARY(Gpufuncs ${LIB_TYPE}
  ${PROJECT_SOURCE_DIR}/Gpufuncs/GpuNegate.cu
  ${PROJECT_SOURCE_DIR}/Gpufuncs/gramsch_update_psi_gpu.cu
  ${PROJECT_SOURCE_DIR}/Gpufuncs/app8_del2_gpu.cu
  ${PROJECT_SOURCE_DIR}/Gpufuncs/GpuEleMul.cu
  ${PROJECT_SOURCE_DIR}/Gpufuncs/GpuFill.cu)

endif(RMG_CUDA_ENABLED)

if(RMG_HIP_ENABLED)
set(MY_HIPCC_OPTIONS -fPIE)
set_source_files_properties(GpuEleMul.hip.cpp PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
hip_add_library(Gpufuncs GpuEleMul.hip.cpp HIPCC_OPTIONS ${MY_HIPCC_OPTIONS})
include_directories("${RMG_SOURCE_DIR}/Headers/")
endif(RMG_HIP_ENABLED)

