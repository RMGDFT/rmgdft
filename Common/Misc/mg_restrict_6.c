/************************** SVN Revision Information **************************
 **    $Id: mg_restrict_6.c 1168 2010-12-03 23:53:46Z btan $    **
******************************************************************************/

/****f* QMD-MGDFT/mg_restrict_6.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void mg_restrict(REAL *full, REAL *half, int dimx, int dimy, int dimz)
 *   get half dimensioned array in the corse grid from full dimensioned 
 *   array in the fine grid. The returned values are smoothed.
 * INPUTS
 *   full[(dimx+10)*(dimy+10)*(dimz+10)]: array in the fine grid
 *      image value must be there
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * OUTPUT
 *   half[(dimx/2)*(dimy/2)*(dimz/2)] array in the corse grid
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>





void mg_restrict_6 (REAL * full, REAL * half, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incz, incy, incx, incz2, incy2, incx2, incx3;
    REAL scale, a0, a1, a2, a3, a4, a5;
    REAL * fulla;
    REAL * fullb;


    incz = 1;
    incy = dimz + 10;
    incx = (dimz + 10) * (dimy + 10);

    incz2 = 1;
    incy2 = dimz / 2 ;
    incx2 = (dimz / 2) * (dimy / 2);
 
    incx3 = (dimz + 10) * (dimy / 2);   

    a0 = 3.0/7.0;
    a1 = 3.0/7.0;
    a2 = (-6.0)/35.0;
    a3 = 1.0/35.0;
    a4 = 0.0;
    a5 = 0.0;

   my_malloc (fulla,(FPX0_GRID / 2)*(FPY0_GRID + 10)*(FPZ0_GRID + 10),REAL);
   my_malloc (fullb,(FPX0_GRID / 2)*(FPY0_GRID / 2)*(FPZ0_GRID + 10),REAL);


        for (ix = 0; ix < dimx / 2; ix++)
        {

            for (iy = 0; iy < dimy + 10; iy++)
            {


                for (iz = 0; iz < dimz + 10; iz++)
                {


                    fulla[ix * incx + iy * incy + iz] =
                               a5* full[(2*ix+0) * incx + iy * incy + iz] +
                               a4* full[(2*ix+1) * incx + iy * incy + iz] +
                               a3* full[(2*ix+2) * incx + iy * incy + iz] +
                               a2* full[(2*ix+3) * incx + iy * incy + iz] + 
                               a1* full[(2*ix+4) * incx + iy * incy + iz] +
                               a0* full[(2*ix+5) * incx + iy * incy + iz] +
                               a1* full[(2*ix+6) * incx + iy * incy + iz] +
                               a2* full[(2*ix+7) * incx + iy * incy + iz] +
                               a3* full[(2*ix+8) * incx + iy * incy + iz] +
                               a4* full[(2*ix+9) * incx + iy * incy + iz] +
                               a5* full[(2*ix+10) * incx + iy * incy + iz] ;


                }               /* end for */

            }                   /* end for */

        }                       /* end for */



        for (ix = 0; ix < dimx / 2; ix++)
        {
        
            for (iy = 0; iy < dimy / 2; iy++)
            {
    

                for (iz = 0; iz < dimz + 10; iz++)
                {


                    fullb[ix * incx3 + iy * incy + iz] =
                               a5* fulla[ix * incx + (2*iy+0) * incy + iz] +
                               a4* fulla[ix * incx + (2*iy+1) * incy + iz] +
                               a3* fulla[ix * incx + (2*iy+2) * incy + iz] +
                               a2* fulla[ix * incx + (2*iy+3) * incy + iz] +
                               a1* fulla[ix * incx + (2*iy+4) * incy + iz] +
                               a0* fulla[ix * incx + (2*iy+5) * incy + iz] +
                               a1* fulla[ix * incx + (2*iy+6) * incy + iz] +
                               a2* fulla[ix * incx + (2*iy+7) * incy + iz] +
                               a3* fulla[ix * incx + (2*iy+8) * incy + iz] +
                               a4* fulla[ix * incx + (2*iy+9) * incy + iz] +
                               a5* fulla[ix * incx + (2*iy+10) * incy + iz] ;


                }               /* end for */

            }                   /* end for */

        }                       /* end for */



        for (ix = 0; ix < dimx / 2; ix++)
        {
        
            for (iy = 0; iy < dimy / 2; iy++)
            {
    

                for (iz = 0; iz < dimz / 2; iz++)
                {


                    half[ix * incx2 + iy * incy2 + iz] =
                               a5* fullb[ix * incx3 + iy* incy + 2*iz+0 ] +
                               a4* fullb[ix * incx3 + iy* incy + 2*iz+1 ] +
                               a3* fullb[ix * incx3 + iy* incy + 2*iz+2 ] +
                               a2* fullb[ix * incx3 + iy* incy + 2*iz+3 ] +
                               a1* fullb[ix * incx3 + iy* incy + 2*iz+4 ] +
                               a0* fullb[ix * incx3 + iy* incy + 2*iz+5 ] +
                               a1* fullb[ix * incx3 + iy* incy + 2*iz+6 ] +
                               a2* fullb[ix * incx3 + iy* incy + 2*iz+7 ] +
                               a3* fullb[ix * incx3 + iy* incy + 2*iz+8 ] +
                               a4* fullb[ix * incx3 + iy* incy + 2*iz+9 ] +
                               a5* fullb[ix * incx3 + iy* incy + 2*iz+10 ] ;


                }               /* end for */

            }                   /* end for */

        }                       /* end for */


        my_free (fulla);
        my_free (fullb);




}                               /* end mg_restrict */


/******/
