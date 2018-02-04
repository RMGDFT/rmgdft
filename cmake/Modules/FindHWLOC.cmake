# - Find libhwloc
#
#  HWLOC_LIBRARIES   - List of libraries when using hwloc
#  HWLOC_FOUND       - True if HWLOC found.

if (HWLOC_LIBRARIES)
  # Already in cache, be silent
  set (HWLOC_FIND_QUIETLY TRUE)
endif (HWLOC_LIBRARIES)

find_library (HWLOC_LIBRARIES NAMES hwloc)

find_path (HWLOC_INCLUDES hwloc.h)

# handle the QUIETLY and REQUIRED arguments and set to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (HWLOC DEFAULT_MSG HWLOC_LIBRARIES HWLOC_INCLUDES)

mark_as_advanced (HWLOC_LIBRARIES HWLOC_INCLUDES)


