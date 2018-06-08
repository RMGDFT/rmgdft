# - Find Scalapack
#
#  SCALAPACK_LIBRARIES   - List of libraries when using scalapack.
#  SCALAPACK_FOUND       - True if SCALAPACK found.

if (SCALAPACK_LIBRARIES)
  # Already in cache, be silent
  set (SCALAPACK_FIND_QUIETLY TRUE)
endif (SCALAPACK_LIBRARIES)

find_library (SCALAPACK_LIBRARIES NAMES scalapack)

# if not found try using the MPI library path as a hint
# since it's often included there
if(NOT SCALAPACK_LIBRARIES)
    get_filename_component(RMG_MPI_LIB_PATH ${MPI_C_LIBRARIES} PATH CACHE)
endif()
find_library (SCALAPACK_LIBRARIES NAMES scalapack
    HINTS "${RMG_MPI_LIB_PATH}" "$ENV{SCALAPACK_LIB}"
)

# handle the QUIETLY and REQUIRED arguments and set SCALAPACK_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SCALAPACK DEFAULT_MSG SCALAPACK_LIBRARIES)

mark_as_advanced (SCALAPACK_LIBRARIES)

