/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>


/*This calulates gradient of beta  using finite differences on the coarse grid*/
/*This function is actually sligtly modified get_grad*/
void partial_beta_fdiff (fftw_complex * beptr, int nldim, rmg_double_t * beta_x, rmg_double_t * beta_y,
                         rmg_double_t * beta_z)
{




    int iz, ix, iy, index, index2;
    rmg_double_t t1x, t2x, t1y, t2y, t1z, t2z;
    rmg_double_t *rptr;
    rmg_double_t *wxr, *wyr, *wzr;
    /*CSS0_GRID *SSSptr; */
    int ixs, iys;
    int ix1, iy1;




    ixs = (nldim + 4) * (nldim + 4);
    iys = (nldim + 4);
    ix1 = nldim * nldim;
    iy1 = nldim;
    wxr = beta_x;
    wyr = beta_y;
    wzr = beta_z;

    my_malloc (rptr, (nldim + 4) * (nldim + 4) * (nldim + 4), rmg_double_t);
    /*SSSptr = (CSS0_GRID *)rptr; */


    /*time1 = my_crtc(); */
    /*trade_images2c(f, SSSptr, b_array); */

    /*Pack data into rptr */
    for (ix = 0; ix < nldim + 4; ix++)
    {

        for (iy = 0; iy < nldim + 4; iy++)
        {

            for (iz = 0; iz < nldim + 4; iz++)
            {

                index = ix * ixs + iy * iys + iz;

                if ((ix < 2) || (iy < 2) || (iz < 2) ||
                    (ix > nldim + 1) || (iy > nldim + 1) || (iz > nldim + 1))
                    rptr[index] = ZERO;

                else
                {

                    index2 = (ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2;
                    rptr[index] = beptr[index2].re;

                }               /*end else */

            }                   /*end for(iz = 0;iz < nldim+4;iz++) */
        }                       /*end for(iy = 0;iy < nldim+4;iy++) */
    }                           /*for(ix = 0;ix < nldim+4;ix++) */


    switch (ct.ibrav)
    {
    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:


        t1x = 8.0 / (12.0 * ct.hxgrid * ct.xside);
        t2x = 1.0 / (12.0 * ct.hxgrid * ct.xside);

        t1y = 8.0 / (12.0 * ct.hygrid * ct.yside);
        t2y = 1.0 / (12.0 * ct.hygrid * ct.yside);

        t1z = 8.0 / (12.0 * ct.hzgrid * ct.zside);
        t2z = 1.0 / (12.0 * ct.hzgrid * ct.zside);

        for (ix = 2; ix < nldim + 2; ix++)
        {

            for (iy = 2; iy < nldim + 2; iy++)
            {

                for (iz = 2; iz < nldim + 2; iz++)
                {

                    wxr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2x * rptr[(ix - 2) * ixs + iy * iys + iz] +
                        -t1x * rptr[(ix - 1) * ixs + iy * iys + iz] +
                        t1x * rptr[(ix + 1) * ixs + iy * iys + iz] +
                        -t2x * rptr[(ix + 2) * ixs + iy * iys + iz];

                    wyr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2y * rptr[ix * ixs + (iy - 2) * iys + iz] +
                        -t1y * rptr[ix * ixs + (iy - 1) * iys + iz] +
                        t1y * rptr[ix * ixs + (iy + 1) * iys + iz] +
                        -t2y * rptr[ix * ixs + (iy + 2) * iys + iz];

                    wzr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2z * rptr[ix * ixs + iy * iys + iz - 2] +
                        -t1z * rptr[ix * ixs + iy * iys + iz - 1] +
                        t1z * rptr[ix * ixs + iy * iys + iz + 1] +
                        -t2z * rptr[ix * ixs + iy * iys + iz + 2];



                }               /* end for */
            }                   /* end for */
        }                       /* end for */

        break;

    case HEXAGONAL:

        error_handler ("HEXAGONAL is not implemented");

        /*The problem is that that app_gradc (and this is just hacked version of that function)
         * uses SSSptr and wx, wy and wz which have 3D indices and I am not sure how to remap the
         * indices. The thing is app_gradc uses CSS0_GRID and P0_GRID types which are defined on the 
         * whole grid but I have only nldim*nldim*nldim size, so I would have to define new data types
         * So, if you want to use HEXAGONAL here, it is fairly easy.
         * If you do not know how to change the code, you should probably compile with FDIFF_BETA set to 0*/
        break;


    default:
        error_handler ("Lattice type not implemented");
    }                           /* end switch */




    my_free (rptr);

    /*time2 = my_crtc();
       rmg_timings(APPGRAD_TIME, (time2 - time1)); */



}                               /* end app_grad */


/******/
