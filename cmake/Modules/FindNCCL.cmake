# - Find NCCL
#
#  NCCL_LIBRARIES   - List of libraries when using FFTW.
#  NCCL_FOUND       - True if FFTW found.

if (NCCL_LIBRARIES)
  # Already in cache, be silent
  set (NCCL_FIND_QUIETLY TRUE)
endif (NCCL_LIBRARIES)

find_library (NCCL_LIBRARIES NAMES nccl)
find_path (NCCL_INCLUDES nccl.h)

# handle the QUIETLY and REQUIRED arguments and set NCCL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NCCL DEFAULT_MSG NCCL_LIBRARIES)

mark_as_advanced (NCCL_LIBRARIES)

