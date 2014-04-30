add_library (RMG_Spin_XC
xcgga.c
xcgga_spin.c
xclda.c                                                                                                                                         
xclsda_spin.c
gcxcblyp.c
gcxcpw91.c
gcxcblyp_spin.c
gcxcpbe.c
lsd_glyp.c
lyp.c
gcxcpw91_spin.c
pbex.c
perdew86.c
gcxcpbe_spin.c
gcxbcp.c
gcxbcp_spin.c
xclda_pz81.c
)
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
include_directories("${RMG_BINARY_DIR}/RMG/Headers/")

