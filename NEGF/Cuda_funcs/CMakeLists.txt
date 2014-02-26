CUDA_ADD_LIBRARY (NEGF_Cudafunc
matrix_inverse_anyprobe_cuda.c
matrix_inverse_cuda.c
rho_munu_cuda.c
Sgreen_cuda.c
genvpsi_gpu.cu
rho_psi_times_psi.cu
)
include_directories("${RMG_SOURCE_DIR}/NEGF/Headers/")
include_directories("${RMG_BINARY_DIR}/NEGF/Headers/")
