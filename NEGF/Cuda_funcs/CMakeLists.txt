ADD_LIBRARY (NEGF_Cudafunc
rho_munu_cuda.c
Sgreen_cuda.c
)
include_directories("${RMG_SOURCE_DIR}/NEGF/Headers/")
include_directories("${RMG_BINARY_DIR}/NEGF/Headers/")
