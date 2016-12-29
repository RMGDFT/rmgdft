# - Find libnuma
#
#  NUMA_LIBRARIES   - List of libraries when using numa
#  NUMA_FOUND       - True if NUMA found.

if (NUMA_LIBRARIES)
  # Already in cache, be silent
  set (NUMA_FIND_QUIETLY TRUE)
endif (NUMA_LIBRARIES)

find_library (NUMA_LIBRARIES NAMES numa)

find_path (NUMA_INCLUDES numa.h)

# handle the QUIETLY and REQUIRED arguments and set to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NUMA DEFAULT_MSG NUMA_LIBRARIES NUMA_INCLUDES)

mark_as_advanced (NUMA_LIBRARIES NUMA_INCLUDES)


