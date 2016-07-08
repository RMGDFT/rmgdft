# - Find Libxc
#
#  LIBXC_INCLUDES    - where to find xc.h
#  LIBXC_LIBRARIES   - List of libraries when using LIBXC.
#  LIBXC_FOUND       - True if LIBXC found.

if (LIBXC_INCLUDES)
  # Already in cache, be silent
  set (LIBXC_FIND_QUIETLY TRUE)
endif (LIBXC_INCLUDES)

find_path (LIBXC_INCLUDES xc.h 
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/include")

find_library (LIBXC_LIBRARIES
NAMES xc
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib64"
"${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib"
)

if(NOT LIBXC_LIBRARIES)
    find_library (LIBXC_LIBRARIES NAMES libxc.a
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib64"
"${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib"
)
endif(NOT LIBXC_LIBRARIES)


find_library (LIBXC_LIBRARIES_F90
NAMES xcf90
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib64"
"${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib"
)

if(NOT LIBXC_LIBRARIES_F90)
    find_library (LIBXC_LIBRARIES_F90 NAMES libxcf90.a
HINTS "${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib64"
"${PROJECT_SOURCE_DIR}/lib/libxc-2.2.2/lib64/lib"
)
endif(NOT LIBXC_LIBRARIES_F90)

# handle the QUIETLY and REQUIRED arguments and set LIBXC_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LIBXC DEFAULT_MSG LIBXC_LIBRARIES LIBXC_INCLUDES)

mark_as_advanced (LIBXC_LIBRARIES LIBXC_INCLUDES)

