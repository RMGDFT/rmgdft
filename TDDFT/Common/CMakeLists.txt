add_library (tddft_Common

read_rhomatrix.cpp
dipole_calculation.cpp
dot_product_orbit_xyz_orbit.cpp
orbit_xyz_orbit.cpp
get_phi_xyz_phi.cpp
mat_global_to_dist.cpp
mat_dist_to_global.cpp
update_TDDFT.cpp
init_TDDFT.cpp
)

include_directories("${RMG_SOURCE_DIR}/TDDFT/Headers/")
include_directories("${RMG_BINARY_DIR}/TDDFT/Headers/")

