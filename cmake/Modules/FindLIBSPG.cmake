# - Find Libspg
#
#  LIBSPG_INCLUDES    - where to find fftw.h
#  LIBSPG_LIBRARIES   - List of libraries when using LIBSPG.
#  LIBSPG_FOUND       - True if LIBSPG found.

if (LIBSPG_LIBRARIES)
  # Already in cache, be silent
  set (LIBSPG_FIND_QUIETLY TRUE)
endif (LIBSPG_LIBRARIES)

find_library (LIBSPG_LIBRARIES 
NAMES symspg
HINTS "${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64/")

if(NOT LIBSPG_LIBRARIES)
    find_library (LIBSPG_LIBRARIES NAMES libsymspg.a
HINTS "${PROJECT_SOURCE_DIR}/lib/spglib-1.6.3/lib64/")
endif(NOT LIBSPG_LIBRARIES)


# handle the QUIETLY and REQUIRED arguments and set LIBSPG_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LIBSPG DEFAULT_MSG LIBSPG_LIBRARIES)

mark_as_advanced (LIBSPG_LIBRARIES)

