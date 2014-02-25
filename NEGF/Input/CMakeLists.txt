add_library (NEGF_Input
read_atoms.c
read_control.c
read_pdb.c
)
include_directories("${RMG_SOURCE_DIR}/NEGF/Headers/")
include_directories("${RMG_BINARY_DIR}/NEGF/Headers/")
