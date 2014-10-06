# - Find Libxc
#
#  LIBXC_INCLUDES    - where to find fftw.h
#  LIBXC_LIBRARIES   - List of libraries when using LIBXC.
#  LIBXC_FOUND       - True if LIBXC found.

if (LIBXC_INCLUDES)
  # Already in cache, be silent
  set (LIBXC_FIND_QUIETLY TRUE)
endif (LIBXC_INCLUDES)

find_path (LIBXC_INCLUDES xc.h 
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/src")

find_library (LIBXC_LIBRARIES 
NAMES xc
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.0.3/lib64/lib64")

# handle the QUIETLY and REQUIRED arguments and set LIBXC_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LIBXC DEFAULT_MSG LIBXC_LIBRARIES LIBXC_INCLUDES)

mark_as_advanced (LIBXC_LIBRARIES LIBXC_INCLUDES)

