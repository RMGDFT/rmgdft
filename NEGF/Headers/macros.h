/************************** SVN Revision Information **************************
 **    $Id: macros.h 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <stdlib.h>
#include <unistd.h>

/* macros */

#define min(a,b) (((a)>(b)) ? (b) : (a))
#define max(a,b) (((a)>(b)) ? (a) : (b))


#ifndef ERROR_HANDLER_C
#  define error_handler(_MSG_) error_handler(pct.thispe, __FILE__, __LINE__, _MSG_),0
#endif

#define where() printf("\nPE: %d file %s line %d\n\n", pct.thispe, __FILE__, __LINE__), fflush(NULL)


#define progress_tag() printf("[ %3d %3d %4d %8.0f ] %s: ", \
                              ct.md_steps, ct.scf_steps, ct.total_scf_steps, \
                              my_crtc() - ct.time0, __FILE__ )




#define my_fopen(_fhandle_, _filename_, _mode_)\
  ( (pct.thispe == 0) ? printf("\nfopening file '%s' mode '%s'\n", _filename_, _mode_) : 0, \
    (_fhandle_ = fopen(_filename_, _mode_)) == NULL ? \
      printf("**** fopening file '%s':", _filename_), error_handler("can't fopen file") : 0 )

#define my_open(_fhandle_, _filename_, _flags_, _mode_)\
  ( (pct.thispe == 0) ? printf("\nopening file '%s' flags '%s' mode '%s'\n", _filename_, #_flags_, #_mode_) : 0, \
    (_fhandle_ = open(_filename_, _flags_, _mode_)) < 0 ? \
      printf("**** opening file '%s':", _filename_), error_handler("can't open file") : 0 )



