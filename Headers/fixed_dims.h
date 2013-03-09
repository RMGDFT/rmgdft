/* Only set these if you know what you are doing. The coarse grid
   settings in the input file will have no effect if these
   are set and the resulting executable will only work as expected
   if everything else matches. Setting these at compile time instead of
   runtime will make the top level finite difference routines run faster.

   You can set FD_XSIZE,FD_YSIZE and FD_ZSIZE in either make_conf.h or
   by a Makefile definition.
*/

#ifdef FD_XSIZE
  #define         FIXED_XDIM      FD_XSIZE
#else
  #define         FIXED_XDIM      pct.PX0_GRID
#endif
#ifdef FD_YSIZE
  #define         FIXED_YDIM      FD_YSIZE
#else
  #define         FIXED_YDIM      pct.PY0_GRID
#endif
#ifdef FD_ZSIZE
  #define         FIXED_ZDIM      FD_ZSIZE
#else
  #define         FIXED_ZDIM      pct.PZ0_GRID
#endif


