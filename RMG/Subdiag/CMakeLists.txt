add_library (RMG_Subdiag
Subdiag.cpp
Subdiag_Elpa.cpp
Subdiag_Lapack.cpp
Subdiag_Cusolver.cpp
Subdiag_Scalapack.cpp
FoldedSpectrum.cpp
FoldedSpectrumScalapack.cpp
FoldedSpectrumIterator.cpp
FoldedSpectrumSetup.cpp
FoldedSpectrumGSE.cpp
FoldedSpectrumScalapackGSE.cpp
FoldedSpectrumOrtho.cpp
FoldedSpectrumScalapackOrtho.cpp
)
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
include_directories("${RMG_BINARY_DIR}/RMG/Headers/")
include_directories(${MAGMA_INCLUDES})
include_directories(${MAGMA_INCLUDES}/magma)

