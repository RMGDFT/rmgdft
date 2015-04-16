# - Find BLACS
#
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (Z_LIBRARIES)
  # Already in cache, be silent
  set (Z_FIND_QUIETLY TRUE)
endif (Z_LIBRARIES)

find_library (Z_LIBRARIES NAMES z)

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Z DEFAULT_MSG Z_LIBRARIES)

mark_as_advanced (Z_LIBRARIES)

