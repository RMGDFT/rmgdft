# - Find PFFT
# Find the native PFFT3 includes and library
#
#  PFFT_INCLUDES    - where to find pfft.h
#  PFFT_LIBRARIES   - List of libraries when using PFFT3.
#  PFFT_FOUND       - True if PFFT3 found.

if (PFFT_INCLUDES)
  # Already in cache, be silent
  set (PFFT_FIND_QUIETLY TRUE)
endif (PFFT_INCLUDES)

find_path (PFFT_INCLUDES pfft.h
HINTS "/usr/local/include" 
      "${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64/include")


find_library (PFFT_LIBRARIES NAMES libpfft.a
                  HINTS "${PROJECT_SOURCE_DIR}/lib/pfft-master/lib64/lib/"
                        "/usr/local/lib64/"
)

# handle the QUIETLY and REQUIRED arguments and set PFFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PFFT DEFAULT_MSG PFFT_LIBRARIES PFFT_INCLUDES)

mark_as_advanced (PFFT_LIBRARIES PFFT_INCLUDES)

