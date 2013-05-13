/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void mg_prolong_6(rmg_double_t *full, rmg_double_t *half, int dimx, int dimy, int dimz)
 *   get full dimensioned array in the fine grid from half dimensioned 
 *   array in the fine grid. Make sure to do the tradeimagex before prolongation.
 * OUTPUTS
 *   full[(dimx)*(dimy)*(dimz)]: array in the fine grid
 *      image value must be there
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * INPUT
 *   half[(dimx/2+10)*(dimy/2+10)*(dimz/2+10)] array in the corse grid
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>





void mg_prolong_6 (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incz, incy, incx, incz2, incy2, incx2, incx3;
    rmg_double_t scale,  a1, a2, a3, a4, a5;
    rmg_double_t * fulla;
    rmg_double_t * fullb;


    incz = 1;
    incy = dimz / 2 + 10;
    incx = (dimz / 2 + 10) * (dimy / 2 + 10);

    incz2 = 1;
    incy2 = dimz;
    incx2 = dimz  * dimy;
 
    incx3 = (dimz / 2+ 10) * dimy;   

    a1 = 150.0/256.0;
    a2 = (-25.0)/256.0;
    a3 = 3.0/256.0;
    a4 = 0.0;
    a5 = 0.0;

   my_malloc (fulla,(pct.FPX0_GRID)*(pct.FPY0_GRID / 2+ 10)*(pct.FPZ0_GRID / 2 + 10),rmg_double_t);
   my_malloc (fullb,(pct.FPX0_GRID)*(pct.FPY0_GRID)*(pct.FPZ0_GRID / 2 + 10),rmg_double_t);


        for (ix = 0; ix < dimx / 2; ix++)
        {

            for (iy = 0; iy < dimy / 2 + 10; iy++)
            {


                for (iz = 0; iz < dimz / 2 + 10; iz++)
                {


                    fulla[((2 * ix)+1) * incx + iy * incy + iz] =
                               a5* half[(ix+1) * incx + iy * incy + iz] +
                               a4* half[(ix+2) * incx + iy * incy + iz] +
                               a3* half[(ix+3) * incx + iy * incy + iz] +
                               a2* half[(ix+4) * incx + iy * incy + iz] + 
                               a1* half[(ix+5) * incx + iy * incy + iz] +
                               a1* half[(ix+6) * incx + iy * incy + iz] +
                               a2* half[(ix+7) * incx + iy * incy + iz] +
                               a3* half[(ix+8) * incx + iy * incy + iz] +
                               a4* half[(ix+9) * incx + iy * incy + iz] +
                               a5* half[(ix+10) * incx + iy * incy + iz] ;

                    fulla[2 * ix * incx + iy * incy + iz] =
                               half[(ix+5) * incx + iy * incy + iz];

                }               /* end for */

            }                   /* end for */

        }                       /* end for */



        for (ix = 0; ix < dimx ; ix++)
        {
        
            for (iy = 0; iy < dimy / 2; iy++)
            {
    

                for (iz = 0; iz < dimz / 2 + 10; iz++)
                {


                    fullb[ix * incx3 + (2 * iy + 1) * incy + iz] =
                               a5* fulla[ix * incx + (iy+1) * incy + iz] +
                               a4* fulla[ix * incx + (iy+2) * incy + iz] +
                               a3* fulla[ix * incx + (iy+3) * incy + iz] +
                               a2* fulla[ix * incx + (iy+4) * incy + iz] +
                               a1* fulla[ix * incx + (iy+5) * incy + iz] +
                               a1* fulla[ix * incx + (iy+6) * incy + iz] +
                               a2* fulla[ix * incx + (iy+7) * incy + iz] +
                               a3* fulla[ix * incx + (iy+8) * incy + iz] +
                               a4* fulla[ix * incx + (iy+9) * incy + iz] +
                               a5* fulla[ix * incx + (iy+10) * incy + iz] ;

                    fullb[ix * incx3 + 2 * iy * incy + iz] =
                               fulla[ix * incx + (iy+5) * incy + iz];


                }               /* end for */

            }                   /* end for */

        }                       /* end for */



        for (ix = 0; ix < dimx ; ix++)
        {
        
            for (iy = 0; iy < dimy ; iy++)
            {
    

                for (iz = 0; iz < dimz / 2; iz++)
                {


                    full[ix * incx2 + iy * incy2 + 2 * iz + 1] =
                               a5* fullb[ix * incx3 + iy* incy + iz+1 ] +
                               a4* fullb[ix * incx3 + iy* incy + iz+2 ] +
                               a3* fullb[ix * incx3 + iy* incy + iz+3 ] +
                               a2* fullb[ix * incx3 + iy* incy + iz+4 ] +
                               a1* fullb[ix * incx3 + iy* incy + iz+5 ] +
                               a1* fullb[ix * incx3 + iy* incy + iz+6 ] +
                               a2* fullb[ix * incx3 + iy* incy + iz+7 ] +
                               a3* fullb[ix * incx3 + iy* incy + iz+8 ] +
                               a4* fullb[ix * incx3 + iy* incy + iz+9 ] +
                               a5* fullb[ix * incx3 + iy* incy + iz+10 ] ;

                    full[ix * incx2 + iy * incy2 + 2 * iz] =
                               fullb[ix * incx3 + iy* incy + iz+5 ];
                    
                }               /* end for */

            }                   /* end for */

        }                       /* end for */


        my_free (fulla);
        my_free (fullb);




}                               /* end mg_prolong_6 */


/******/
