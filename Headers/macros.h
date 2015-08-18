/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#ifndef MACROS_H_INCLUDED
#define MACROS_H_INCLUDED

#if !(defined(_WIN32) || defined(_WIN64))
    #include <unistd.h>
#else
    #ifndef NOMINMAX
      #define NOMINMAX
    #endif
    #include <io.h>
    #include <windows.h>
    #define sleep(x) Sleep(1000.0*x)

#endif
#include <stdlib.h>

/* macros */

/* Some useful macros */
/*
 * determines position of given pe
 */
#define PE2XYZ(pe,x,y,z) x = pe; z = x % PE_Z; x /= PE_Z; y = x % PE_Y; x /= PE_Y

/*
 * determines pe at given position
 */
#define XYZ2PE(x,y,z,pe) pe = (x)*get_PE_Y()*get_PE_Z() + (y)*get_PE_Z() + (z)

#define DEBUG 0

#if DEBUG
  #if (defined(_WIN32) || defined(_WIN64))
    #define Dprintf( format, args,...) fprintf (stderr, "\n#DEBUG from PE %d in %s line %d:    \t" format "\n", pct.gridpe, __FILE__, __LINE__, ##args), fflush(NULL)
  #else
    #define Dprintf( format, args...) fprintf (stderr, "\n#DEBUG from PE %d in %s line %d:    \t" format "\n", pct.gridpe, __FILE__, __LINE__, ##args), fflush(NULL)
  #endif
#else
  #if (defined(_WIN32) || defined(_WIN64))
    #define Dprintf( format, args,...) ;
  #else
    #define Dprintf( format, args...) ;
  #endif
#endif

//#define dprintf( format, args...) fprintf (stderr, "\ngrid rank %d of spin %d:    \t"format"\n", pct.gridpe, pct.spinpe,  ##args), fflush(NULL), fsync ( fileno (ct.logfile))
#if (defined(_WIN32) || defined(_WIN64))
  #define dprintf( format, args,...) fprintf (stderr, "\nIMG %d/%d:PE %d, GRID %d/%d:PE %d,\t" format "\n", pct.thisimg+1, pct.images, pct.imgpe, pct.spinpe+1, pct.grids, pct.gridpe,  ##args), fflush(NULL)
#else
  #define dprintf( format, args...) fprintf (stderr, "\nIMG %d/%d:PE %d, GRID %d/%d:PE %d,\t" format "\n", pct.thisimg+1, pct.images, pct.imgpe, pct.spinpe+1, pct.grids, pct.gridpe,  ##args), fflush(NULL)
#endif
//#define dprintf( format, args...) fprintf (stderr, "\n#WARNING from IMG PE %d in IMG %d  of grid rank %d of spin %d:    \t"format"\n", pct.imgpe, pct.thisimg, pct.gridpe, pct.thisspin,  ##args), fflush(NULL)


#define my_strncpy(buf1, buf2, num) strncpy(buf1, buf2, num), buf1[num]=0

#define rmg_min(a,b) (((a)>(b)) ? (b) : (a))
#define rmg_max(a,b) (((a)>(b)) ? (a) : (b))

#define even(a)  ( (a) % 2 == 0 )
#define odd(a)   ( (a) % 2 == 1 )


#if (defined(_WIN32) || defined(_WIN64))
    #define printf( message,... ) \
	 ((pct.imgpe == 0) ? fprintf( ct.logfile, message, __VA_ARGS__ ): 0)
//	 ((pct.imgpe == 0) ? fprintf( ct.logfile, message ), fflush (NULL), fsync ( fileno (ct.logfile) ): 0)
#else
    #define printf( message... ) \
	 ((pct.imgpe == 0) ? fprintf( ct.logfile, message ): 0)
//	 ((pct.imgpe == 0) ? fprintf( ct.logfile, message ), fflush (NULL), fsync ( fileno (ct.logfile) ): 0)
#endif
	

/* variadic error_handler, use is the same as printf. Calls MPI_Abort, since MPI_Finalize hangs if called on single PE. */
#if (defined(_WIN32) || defined(_WIN64))
  #define error_handler( message,... ) \
    fprintf (stderr, "\nExit from PE %d of image %d, in file %s, line %d\nPE %d Error Message is: ", pct.gridpe, pct.thisimg+1, __FILE__, __LINE__, pct.gridpe), \
    printf ("\nExit from PE %d of image %d, in file %s, line %d\nPE %d Error Message is: ", pct.gridpe, pct.thisimg+1, __FILE__, __LINE__, pct.gridpe), \
	fprintf (stderr,  message, __VA_ARGS__ ), \
	fprintf (stderr,  "\n\n" ), \
	printf ( message, __VA_ARGS__ ), \
	printf ( "\n\n" ), \
    fflush (NULL), \
    sleep (2), \
	MPI_Abort( MPI_COMM_WORLD, 0 )
#else
  #define error_handler( message... ) \
    fprintf (stderr, "\nExit from PE %d of image %d, in file %s, line %d\nPE %d Error Message is: ", pct.gridpe, pct.thisimg+1, __FILE__, __LINE__, pct.gridpe), \
    printf ("\nExit from PE %d of image %d, in file %s, line %d\nPE %d Error Message is: ", pct.gridpe, pct.thisimg+1, __FILE__, __LINE__, pct.gridpe), \
	fprintf (stderr,  message ), \
	fprintf (stderr,  "\n\n" ), \
	printf ( message ), \
	printf ( "\n\n" ), \
    fflush (NULL), \
	fsync ( fileno (ct.logfile) ), \
    sleep (2), \
	MPI_Abort( MPI_COMM_WORLD, 0 )
#endif

#define progress_tag() printf("[ %3d %3d %4d %8.0f ] %s: ", \
                              ct.md_steps, ct.scf_steps, ct.total_scf_steps, \
                              my_crtc() - ct.time0, __FILE__ )


#define my_fopen(_fhandle_, _filename_, _mode_) do {\
    _fhandle_ = fopen(_filename_, _mode_);\
    if (_fhandle_ == NULL)\
        error_handler("macro my_fopen: can't fopen file %s", _filename_);\
} while (0)
	  


#define my_open(_fhandle_, _filename_, _flags_, _mode_) do {\
    _fhandle_ = open(_filename_, _flags_, _mode_);\
    if (_fhandle_ < 0)\
        error_handler("macro my_open: can't open file %s", _filename_);\
} while (0)

#endif
