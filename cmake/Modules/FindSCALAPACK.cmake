# - Find Scalapack
#
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (SCALAPACK_LIBRARIES)
  # Already in cache, be silent
  set (SCALAPACK_FIND_QUIETLY TRUE)
endif (SCALAPACK_LIBRARIES)

find_library (SCALAPACK_LIBRARIES NAMES scalapack)

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SCALAPACK DEFAULT_MSG SCALAPACK_LIBRARIES)

mark_as_advanced (SCALAPACK_LIBRARIES)

