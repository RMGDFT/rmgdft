# - Find BLACS
#
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (BZ2_LIBRARIES)
  # Already in cache, be silent
  set (BZ2_FIND_QUIETLY TRUE)
endif (BZ2_LIBRARIES)

find_library (BZ2_LIBRARIES NAMES bz2)

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BZ2 DEFAULT_MSG BZ2_LIBRARIES)

mark_as_advanced (BZ2_LIBRARIES)

