/****f* QMD-MGDFT/mg_prolong_MAX10.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 * FUNCTION
 *   void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int grid_ratio, int order)
 *   get full dimensioned array in the fine grid from half dimensioned 
 *   array in the fine grid. 
 *   Attention: order <= Max_order = 10
 * OUTPUTS
 *   full[(dimx)*(dimy)*(dimz)]: array in the fine grid
 *      image value must be there
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * INPUT
 *   half[(dimx/grid_ratio+10)*(dimy/grid_ratio+10)*(dimz/grid_ratio+10)] array in the corse grid
 * PARENTS
 *   get_rho.c or get_new_rho.c
 * CHILDREN
 *   cgen_prolong.c 
 * SOURCE
 */



#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"
#include <TradeImages.h>


#define MAX_ORDER 10

void mg_prolong_MAX10 (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz, int grid_ratio, int order)
{

    int ix, iy, iz, i, k;
    int incy, incx,  incy2, incx2, incx3;
    double a[MAX_ORDER][MAX_ORDER];
    double *fulla, *sg_half;
    double *fullb;
    double c[MAX_ORDER];
    double fraction;

    sg_half = new double[(half_dimx + 10) * (half_dimy + 10) * (half_dimz + 10)];
    trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, 5, FULL_TRADE);


    incy = dimz / grid_ratio + 10;
    incx = (dimz / grid_ratio + 10) * (dimy / grid_ratio + 10);

    incy2 = dimz;
    incx2 = dimz * dimy;

    incx3 = (dimz / grid_ratio + 10) * dimy;

    /*Order has to be even number */
    if (order % 2)
        error_handler ("This function works only for even orders, but order %d was specified",
                       order);


    for (ix = 0; ix < MAX_ORDER; ix++)
    {
        for (iy = 0; iy < MAX_ORDER; iy++)
        {

            a[ix][iy] = 0.0;
        }
    }


    for (i = 0; i < grid_ratio; i++)
    {

        fraction = (double) i / grid_ratio;
        cgen_prolong (c, fraction, order);

        for (iy = 0; iy < order; iy++)
        {
            k = iy + (MAX_ORDER - order) / 2;
            a[i][k] = c[iy];
            //  printf("  a[%d][%d]= %f \n ", i, k, a[i][k]);
        }




    }

    fulla = new double[dimx * (dimy / grid_ratio + 10) * (dimz / grid_ratio + 10)];
    fullb = new double[ dimx * dimy * (dimz / grid_ratio + 10)];


    for (ix = 0; ix < dimx / grid_ratio; ix++)
    {

        for (iy = 0; iy < dimy / grid_ratio + 10; iy++)
        {


            for (iz = 0; iz < dimz / grid_ratio + 10; iz++)
            {


                for (i = 0; i < grid_ratio; i++)


                {


                    fulla[((grid_ratio * ix) + i) * incx + iy * incy + iz] =
                        a[i][0] * sg_half[(ix + 1) * incx + iy * incy + iz] +
                        a[i][1] * sg_half[(ix + 2) * incx + iy * incy + iz] +
                        a[i][2] * sg_half[(ix + 3) * incx + iy * incy + iz] +
                        a[i][3] * sg_half[(ix + 4) * incx + iy * incy + iz] +
                        a[i][4] * sg_half[(ix + 5) * incx + iy * incy + iz] +
                        a[i][5] * sg_half[(ix + 6) * incx + iy * incy + iz] +
                        a[i][6] * sg_half[(ix + 7) * incx + iy * incy + iz] +
                        a[i][7] * sg_half[(ix + 8) * incx + iy * incy + iz] +
                        a[i][8] * sg_half[(ix + 9) * incx + iy * incy + iz] +
                        a[i][9] * sg_half[(ix + 10) * incx + iy * incy + iz];

                }               /* end for */

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy / grid_ratio; iy++)
        {


            for (iz = 0; iz < dimz / grid_ratio + 10; iz++)
            {


                for (i = 0; i < grid_ratio; i++)


                {

                    fullb[ix * incx3 + (grid_ratio * iy + i) * incy + iz] =
                        a[i][0] * fulla[ix * incx + (iy + 1) * incy + iz] +
                        a[i][1] * fulla[ix * incx + (iy + 2) * incy + iz] +
                        a[i][2] * fulla[ix * incx + (iy + 3) * incy + iz] +
                        a[i][3] * fulla[ix * incx + (iy + 4) * incy + iz] +
                        a[i][4] * fulla[ix * incx + (iy + 5) * incy + iz] +
                        a[i][5] * fulla[ix * incx + (iy + 6) * incy + iz] +
                        a[i][6] * fulla[ix * incx + (iy + 7) * incy + iz] +
                        a[i][7] * fulla[ix * incx + (iy + 8) * incy + iz] +
                        a[i][8] * fulla[ix * incx + (iy + 9) * incy + iz] +
                        a[i][9] * fulla[ix * incx + (iy + 10) * incy + iz];

                }               /* end for */

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {


            for (iz = 0; iz < dimz / grid_ratio; iz++)
            {


                for (i = 0; i < grid_ratio; i++)


                {


                    full[ix * incx2 + iy * incy2 + grid_ratio * iz + i] =
                        a[i][0] * fullb[ix * incx3 + iy * incy + iz + 1] +
                        a[i][1] * fullb[ix * incx3 + iy * incy + iz + 2] +
                        a[i][2] * fullb[ix * incx3 + iy * incy + iz + 3] +
                        a[i][3] * fullb[ix * incx3 + iy * incy + iz + 4] +
                        a[i][4] * fullb[ix * incx3 + iy * incy + iz + 5] +
                        a[i][5] * fullb[ix * incx3 + iy * incy + iz + 6] +
                        a[i][6] * fullb[ix * incx3 + iy * incy + iz + 7] +
                        a[i][7] * fullb[ix * incx3 + iy * incy + iz + 8] +
                        a[i][8] * fullb[ix * incx3 + iy * incy + iz + 9] +
                        a[i][9] * fullb[ix * incx3 + iy * incy + iz + 10];


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    delete [] fulla;
    delete [] fullb;
    delete [] sg_half;

}                               /* end mg_prolong_MAX10 */


/******/
