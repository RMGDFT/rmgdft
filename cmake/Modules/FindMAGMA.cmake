# - Find MAGMA
# Find the native MAGMA includes and library
#
#  MAGMA_INCLUDES    - where to find magma.h
#  MAGMA_LIBRARIES   - List of libraries when using MAGMA.
#  MAGMA_FOUND       - True if MAGMA found.

if (MAGMA_INCLUDES)
  # Already in cache, be silent
  set (MAGMA_FIND_QUIETLY TRUE)
endif (MAGMA_INCLUDES)

find_path (MAGMA_INCLUDES magma.h)
if(NOT MAGMA_INCLUDES)
    find_path (MAGMA_INCLUDES magma/magma.h)
endif(NOT MAGMA_INCLUDES)

find_library (MAGMA_LIBRARIES NAMES magma PATH_SUFFIXES magma)
find_library (MAGMABLAS_LIBRARIES NAMES magmablas PATH_SUFFIXES magmablas)

# handle the QUIETLY and REQUIRED arguments and set MAGMA_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (MAGMA DEFAULT_MSG MAGMA_LIBRARIES MAGMA_INCLUDES)
find_package_handle_standard_args (MAGMABLAS DEFAULT_MSG MAGMABLAS_LIBRARIES MAGMA_INCLUDES)
#find_package_handle_standard_args (MAGMA DEFAULT_MSG MAGMA_LIBRARIES
#MAGMA_INCLUDES MAGMABLAS_LIBRARIES)

mark_as_advanced (MAGMA_LIBRARIES MAGMA_INCLUDES)

