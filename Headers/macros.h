#pragma once
#ifndef MACROS_H_INCLUDED
#define MACROS_H_INCLUDED

#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>


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
    #define Dprintf( format, args...) fprintf (stderr, "\n#DEBUG from PE %d in %s line %d:    \t" format "\n", pct.gridpe, __FILE__, __LINE__, ##args), fflush(NULL)
#else
    #define Dprintf( format, args...) ;
#endif

#define dprintf( format, args...) fprintf (stderr, "\nIMG %d/%d:PE %d, GRID %d/%d:PE %d,\t" format "\n", pct.thisimg+1, pct.images, pct.imgpe, pct.spinpe+1, pct.grids, pct.gridpe,  ##args), fflush(NULL)


#define my_fopen(_fhandle_, _filename_, _mode_) do {\
    _fhandle_ = fopen(_filename_, _mode_);\
    if (_fhandle_ == NULL){\
        printf("macro my_fopen: can't fopen file %s", _filename_);\
        fflush (NULL);\
        raise(SIGTERM);\
    }\
} while (0)
	  


#define my_open(_fhandle_, _filename_, _flags_, _mode_) do {\
    _fhandle_ = open(_filename_, _flags_, _mode_);\
    if (_fhandle_ < 0){\
        printf("macro my_open: can't open file %s", _filename_);\
        printf("macro my_fopen: can't fopen file %s", _filename_);\
        fflush (NULL);\
        raise(SIGTERM);\
    }\
} while (0)

#endif
