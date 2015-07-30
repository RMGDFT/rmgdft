# - Find PMRRR
#

if (PMRRR_LIBRARIES)
  # Already in cache, be silent
  set (PMRRR_FIND_QUIETLY TRUE)
endif (PMRRR_LIBRARIES)

find_library (PMRRR_LIBRARIES NAMES pmrrr
HINTS "${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build/external/pmrrr")

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PMRRR DEFAULT_MSG PMRRR_LIBRARIES)

mark_as_advanced (PMRRR_LIBRARIES)


