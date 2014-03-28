# - Find BLACS
#
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (BLACS_LIBRARIES)
  # Already in cache, be silent
  set (BLACS_FIND_QUIETLY TRUE)
endif (BLACS_LIBRARIES)

find_library (BLACS_LIBRARIES NAMES mpiblacs)
find_library (BLACSCINIT_LIBRARIES NAMES mpiblacsCinit)

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BLACS DEFAULT_MSG BLACS_LIBRARIES BLACSCINIT_LIBRARIES ) 

mark_as_advanced (BLACS_LIBRARIES BLACSCINIT_LIBRARIES )

