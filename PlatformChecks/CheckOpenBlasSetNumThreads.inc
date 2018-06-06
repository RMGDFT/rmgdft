set(CMAKE_REQUIRED_INCLUDES "${RMG_BINARY_DIR}/Headers;${RMG_SOURCE_DIR}/Headers")
set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
CHECK_C_SOURCE_COMPILES("int main(int argc, char **argv) {
                           openblas_set_num_threads(2);
                         }" RMG_OPENBLAS_SET_NUM_THREADS_INTERNAL)

