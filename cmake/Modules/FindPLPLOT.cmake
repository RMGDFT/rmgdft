# Find PLplot header and library.

# This module defines the following uncached variables:
#  PLPLOT_FOUND, if false, do not try to use PLplot.
#  PLplot_INCLUDE_DIR, where to find plplot.h.
#  PLplot_LIBRARIES, the libraries to link against to use PLplot
#  PLplot_LIBRARY_DIR, the directory where the PLplot library is found.

if (PLplot_INCLUDE_DIR)
  # Already in cache, be silent
  set (PLplot_FIND_QUIETLY TRUE)
endif (PLplot_INCLUDE_DIR)

find_path (PLplot_INCLUDE_DIR plplot.h)
if(NOT PLplot_INCLUDE_DIR)
    find_path (PLplot_INCLUDE_DIR plplot/plplot.h)
endif(NOT PLplot_INCLUDE_DIR)

find_library (PLplot_LIBRARY NAMES plplot)
find_library (PLplot_cxx_LIBRARY NAMES plplotcxx)

# handle the QUIETLY and REQUIRED arguments and set PLPLOT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PLPLOT DEFAULT_MSG PLplot_LIBRARY PLplot_cxx_LIBRARY PLplot_INCLUDE_DIR)

mark_as_advanced (PLplot_LIBRARY PLplot_cxx_LIBRARY PLplot_INCLUDE_DIR)

