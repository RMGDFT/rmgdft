add_library (US_PP
InitSpinOrbit.cpp
PackGftoc.cpp
GetPhaseSpecies.cpp
InitClebschGordan.cpp
Ylm.cpp
CubicHarmonic.cpp
InitLocalizedWeight.cpp
InitDelocalizedWeight.cpp
InitDelocalizedWeight_onek.cpp
DelocalizedWeight_one.cpp
InitDelocalizedOrbital.cpp
InitLocalizedOrbital.cpp
InitWeightOne.cpp
InitDelocalizedObject.cpp
InitSemilocalBessel.cpp
#PackFine2Rhogrid.cpp
FindFftwPhaseLocalpp.cpp
Radius2grid.cpp
FindPhase.cpp
FindPhaseKpoint.cpp
Atomic.cpp
InitPseudo.cpp
InitLocalObject.cpp
init_efield.cpp
#weight_shift_center.c
ylmr2.cpp
linint.cpp
get_QnmL.cpp
GetQI.cpp
get_qqq.cpp
get_qqq_dk.cpp
radint1.cpp
get_ddd.cpp
radint.cpp
)
include_directories("${RMG_SOURCE_DIR}/Headers/")

