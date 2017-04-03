# - Find FFTW
# Find the native FFTW3 includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW3.
#  FFTW_FOUND       - True if FFTW3 found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)
if(NOT FFTW_INCLUDES)
    find_path (FFTW_INCLUDES dfftw3.h)
endif(NOT FFTW_INCLUDES)

find_library (FFTW_LIBRARIES NAMES dfftw3)
if(NOT FFTW_LIBRARIES)
    find_library (FFTW_LIBRARIES NAMES libfftw3.a)
endif(NOT FFTW_LIBRARIES)

if(NOT FFTW_LIBRARIES)
    find_library (FFTW_LIBRARIES NAMES fftw3)
endif(NOT FFTW_LIBRARIES)

# use libfftw3_mpi.so if you have errors like "fftw_mktensor_4d"
#find_library (FFTW_MPI_LIBRARIES NAMES libfftw3_mpi.so)
#find_library (FFTW_MPI_LIBRARIES NAMES fftw3_mpi)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

