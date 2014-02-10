add_library (MG
mg_restrict.c
mg_sizes.c
mg_restrict_f.c
mgrid_solv.c
mgrid_solv_f.c
mg_prolong_f.c
get_vh.c
mg_prolong.c
find_grid_owner.c
)
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
