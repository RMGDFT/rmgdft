add_library (RMG_Subdiag
ApplyOperators.cpp
Subdiag.cpp
Subdiag_Elpa.cpp
Subdiag_Lapack.cpp
Subdiag_Magma.cpp
Subdiag_Scalapack.cpp
FoldedSpectrum.cpp
FoldedSpectrumScalapack.cpp
FoldedSpectrumIterator.cpp
FoldedSpectrumSetup.cpp
Rmg_dsygvd_gpu.cpp
Rmg_zhegvd_gpu.cpp
FoldedSpectrumGSE.cpp
FoldedSpectrumScalapackGSE.cpp
FoldedSpectrumOrtho.cpp
FoldedSpectrumScalapackOrtho.cpp
)
include_directories("${RMG_SOURCE_DIR}/RMG/Headers/")
include_directories("${RMG_BINARY_DIR}/RMG/Headers/")
include_directories(${MAGMA_INCLUDES})
include_directories(${MAGMA_INCLUDES}/magma)

