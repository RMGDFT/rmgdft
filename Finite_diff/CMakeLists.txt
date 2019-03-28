add_library (Finite_diff
ApplyGradient.cpp
AppGradPfft.cpp
ApplyLaplacian.cpp
ApplyAOperator.cpp
ApplyBOperator.cpp
AppCir.cpp
AppCilrFourth.cpp
AppCilrSixth.cpp
AppCilrDriver.cpp
FiniteDiff_Cbindings.cpp
app_grad.c
app6_del2.c
)
include_directories("${RMG_SOURCE_DIR}/Headers/")

