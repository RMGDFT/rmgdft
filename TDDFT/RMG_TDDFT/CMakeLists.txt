add_library (RmgTddft
RmgTddft.cpp
HmatrixUpdate.cpp
HSmatrix.cpp
GetNewRho_rmgtddft.cpp
)
include_directories("${RMG_SOURCE_DIR}/TDDFT/Headers/")
include_directories("${RMG_BINARY_DIR}/TDDFT/Headers/")

