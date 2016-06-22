add_library (Input
AutoSet.cpp
SetupGrids.cpp
GetPrimeFactors.cpp
ReadForces.cpp
ReadVelocities.cpp
ReadInit.cpp
ReadRmgAtoms.cpp
ReadDynamics.cpp
ReadKpoints.cpp
ReadKpointsBandstructure.cpp
Verify.cpp
InputKey.cpp
CheckValue.cpp
RmgInputFile.cpp
ReadCommon.cpp
MapElements.cpp
Upf.cpp
ReadPseudo.cpp
ReadTFAtoms.cpp
)
include_directories("${RMG_SOURCE_DIR}/Headers/")
if(OPENBABEL_FOUND)
get_filename_component(OPENBABEL_PARENT_DIR ${OPENBABEL_INCLUDES} PATH)
include_directories("${OPENBABEL_PARENT_DIR}")
endif(OPENBABEL_FOUND)

