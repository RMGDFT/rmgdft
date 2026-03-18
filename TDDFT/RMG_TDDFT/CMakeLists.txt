add_library (RmgTddft
rmg_tddft.cpp
init_point_charge_pot.cpp
RmgTddft.cpp
HmatrixUpdate.cpp
HSmatrix.cpp
GetNewRho_rmgtddft.cpp
VecPmatrix.cpp
CurrentNlpp.cpp
CurrentOperator.cpp
EpsilonMatrix.cpp
MatDiagSet.cpp
TddftEnergy.cpp
)
include_directories("${RMG_SOURCE_DIR}/TDDFT/Headers/")
include_directories("${RMG_BINARY_DIR}/TDDFT/Headers/")

