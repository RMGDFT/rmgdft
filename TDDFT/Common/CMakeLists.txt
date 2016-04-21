add_library (tddft_Common

read_rhomatrix.c
dipole_calculation.c
dot_product_orbit_xyz_orbit.c
orbit_xyz_orbit.c
get_phi_xyz_phi.c
mat_global_to_dist.c
mat_dist_to_global.c
update_TDDFT.c
init_TDDFT.c
)

add_library (RmgTddft
WriteData_rmgtddft.cpp
ReadData_rmgtddft.cpp
dipole_calculation.c
mat_global_to_dist.c
mat_dist_to_global.c
RmgTddft.cpp
../../ON/Common/VhDriver.cpp
HmatrixUpdate.cpp
HSmatrix.cpp
GetNewRho_rmgtddft.cpp
)
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_BINARY_DIR}/Headers/")

