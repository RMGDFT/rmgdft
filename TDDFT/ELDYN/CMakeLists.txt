add_library (Eldyn
Magnus.cpp
initmd.f90
timing.f90
utilmod2.f90
commutp.f90
diagev.f90
dochol.f90
eldyn.f90
)
add_dependencies(Eldyn eldyn_mod)
add_library (eldyn_mod
initmd.f90
timing.f90
utilmod2.f90
)

