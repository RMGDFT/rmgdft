# - Find Openbabel
#
#  OPENBABEL_INCLUDES    - where to find fftw3.h
#  OPENBABEL_LIBRARIES   - List of libraries when using OPENBABEL.
#  OPENBABEL_FOUND       - True if OPENBABEL found.

if (OPENBABEL_INCLUDES)
  # Already in cache, be silent
  set (OPENBABEL_FIND_QUIETLY TRUE)
endif (OPENBABEL_INCLUDES)

#find_path (OPENBABEL_INCLUDES babelconfig.h HINTS "/usr/include/openbabel-2.0/openbabel")
find_path (OPENBABEL_INCLUDES babelconfig.h HINTS "/usr/include/openbabel" "/usr/include/openbabel-2.0/openbabel")
find_library (OPENBABEL_LIBRARIES NAMES openbabel)
if(NOT OPENBABEL_LIBRARIES)
    find_library (OPENBABEL_LIBRARIES NAMES libopenbabel.so)
endif(NOT OPENBABEL_LIBRARIES)

# handle the QUIETLY and REQUIRED arguments and set OPENBABEL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OPENBABEL DEFAULT_MSG OPENBABEL_LIBRARIES OPENBABEL_INCLUDES)

mark_as_advanced (OPENBABEL_LIBRARIES OPENBABEL_INCLUDES)

