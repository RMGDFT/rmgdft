add_library(dftd3 STATIC 
dftd3_driver.f90
lib/api.f90
lib/common.f90
lib/core.f90
lib/pars.f90
lib/sizes.f90
RmgDftd3.cpp
)
