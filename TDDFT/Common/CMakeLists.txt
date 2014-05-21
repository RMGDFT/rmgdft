add_library (tddft_Common
dot_product_orbit_xyz_orbit.c
orbit_xyz_orbit.c
get_phi_xyz_phi.c
update_TDDFT.c
init_TDDFT.c
)
include_directories("${RMG_SOURCE_DIR}/ON/Headers/")
include_directories("${RMG_BINARY_DIR}/ON/Headers/")
