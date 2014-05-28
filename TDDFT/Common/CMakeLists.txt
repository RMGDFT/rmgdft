add_library (tddft_Common
RmgPrintTimings.cpp
dipole_calculation.c
dot_product_orbit_xyz_orbit.c
orbit_xyz_orbit.c
get_phi_xyz_phi.c
mat_global_to_dist.c
mat_dist_to_global.c
update_TDDFT.c
init_TDDFT.c
)
include_directories("${RMG_SOURCE_DIR}/Headers/")
include_directories("${RMG_BINARY_DIR}/Headers/")

