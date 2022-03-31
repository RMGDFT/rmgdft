add_library (Diagonalization
DsyevjDriver.cpp
DsyevdDriver.cpp
DsygvdDriver.cpp
DsygvdMgDriver.cpp
DsygvjDriver.cpp
ZhegvdDriver.cpp
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

