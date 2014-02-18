add_library (ON_Input
read_atoms.c
read_control.c
read_pdb.c
)
include_directories("${RMG_SOURCE_DIR}/ON/Headers/")
include_directories("${RMG_BINARY_DIR}/ON/Headers/")
