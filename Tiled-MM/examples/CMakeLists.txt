add_executable(multiply "multiply.cpp")
target_link_libraries(multiply PRIVATE Tiled-MM cxxopts::cxxopts)

if(TILEDMM_GPU_BACKEND STREQUAL "CUDA")
    add_executable(cublasXt-multiply "cublasXt-multiply.cpp")
    target_link_libraries(cublasXt-multiply PRIVATE Tiled-MM cxxopts::cxxopts)
endif()
