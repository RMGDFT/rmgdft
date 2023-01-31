add_executable(multiply "multiply.cpp")
target_link_libraries(multiply Tiled-MM cxxopts)

if(TILEDMM_GPU_BACKEND STREQUAL "CUDA")
    add_executable(cublasXt-multiply "cublasXt-multiply.cpp")
    target_link_libraries(cublasXt-multiply Tiled-MM cxxopts)
endif()
