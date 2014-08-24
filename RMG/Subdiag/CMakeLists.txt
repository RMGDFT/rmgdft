add_library (RMG_Subdiag
ApplyOperators.cpp
Subdiag.cpp
SubdiagGemm.cpp
Subdiag_Lapack.cpp
Subdiag_Magma.cpp
Subdiag_Scalapack.cpp
)
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
include_directories("${RMG_BINARY_DIR}/RMG/Headers/")

