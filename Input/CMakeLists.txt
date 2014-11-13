add_library (Input
ReadForces.cpp
ReadVelocities.cpp
ReadInit.cpp
ReadRmgAtoms.cpp
ReadDynamics.cpp
ReadKpoints.cpp
Verify.cpp
InputKey.cpp
CheckValue.cpp
RmgInputFile.cpp
ReadCommon.cpp
MapElements.cpp
Upf.cpp
ReadPseudo.cpp
)
include_directories("${RMG_SOURCE_DIR}/Headers/")
if(OPENBABEL_FOUND)
get_filename_component(OPENBABEL_PARENT_DIR ${OPENBABEL_INCLUDES} PATH)
include_directories("${OPENBABEL_PARENT_DIR}")
endif(OPENBABEL_FOUND)

