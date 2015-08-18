# - Find ELEMENTAL
#

if (ELEMENTAL_LIBRARIES)
  # Already in cache, be silent
  set (ELEMENTAL_FIND_QUIETLY TRUE)
endif (ELEMENTAL_LIBRARIES)

find_path (ELEMENTAL_INCLUDES El.hpp
HINTS "${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build/include")
find_library (ELEMENTAL_LIBRARIES NAMES El
HINTS "${PROJECT_SOURCE_DIR}/lib/Elemental-0.85/build")

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (ELEMENTAL DEFAULT_MSG
ELEMENTAL_LIBRARIES ELEMENTAL_INCLUDES)

mark_as_advanced (ELEMENTAL_LIBRARIES)


