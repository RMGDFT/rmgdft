/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#include "main.h"
#include "common_prototypes.h"
#include <float.h>
#include <math.h>





void mg_prolong_6 (double * full, double * half, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incz, incy, incx, incz2, incy2, incx2, incx3;
    double scale,  a1, a2, a3, a4, a5;
    double * fulla;
    double * fullb;


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

   my_malloc (fulla,(get_FPX0_GRID())*(get_FPY0_GRID() / 2+ 10)*(get_FPZ0_GRID() / 2 + 10),double);
   my_malloc (fullb,(get_FPX0_GRID())*(get_FPY0_GRID())*(get_FPZ0_GRID() / 2 + 10),double);


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
