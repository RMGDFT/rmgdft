add_library (OrbitalProjection
#get_Hij_update.c
#get_new_rho_local.c
GetNewRhoLocal.cpp
HijUpdateNCpp.cpp
InitStatedistribute.cpp
MapOrbitalToProcess.cpp
MatrixToLocal.cpp
MatrixToGlobal.cpp
)
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_BINARY_DIR}/Headers/")

