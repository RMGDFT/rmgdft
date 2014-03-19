/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "FiniteDiff.h"
#include "TradeImages.h"


rmg_double_t app_cilr_sixth (rmg_double_t * psi, rmg_double_t *a_psi, rmg_double_t *b_psi, rmg_double_t *vtot_eig_s, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int numgrid, tid, P0_BASIS;
    rmg_double_t *rptr;

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    rmg_double_t ecxy, ecxz, ecyz, cc, fcx, fcy, fcz, cor;
    rmg_double_t fc2x, fc2y, fc2z, tcx, tcy, tcz;
    rmg_double_t ihx, ihy, ihz, retv;
    rmg_double_t c000, c100, c110, c200;
    int pbasis = dimx * dimy * dimz, itid;
    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);

    P0_BASIS = get_P0_BASIS();

    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;


    my_malloc (rptr, sbasis + 64, rmg_double_t);
    trade_imagesx (psi, rptr, dimx, dimy, dimz, 2, FULL_TRADE);

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

    cc = (-116.0 / 90.0) * (ihx + ihy + ihz);

    fcx = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihx;
    fcy = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihy;
    fcz = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihz;

    ecxy = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihz;
    ecxz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihy;
    ecyz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihx;

    cor = (1.0 / 144.0) * (ihx + ihy + ihz);

    fc2x = (1.0 / 120.0) * (ihy + ihz);
    fc2y = (1.0 / 120.0) * (ihx + ihz);
    fc2z = (1.0 / 120.0) * (ihx + ihy);

    tcx = (-1.0 / 240.0) * ihx;
    tcy = (-1.0 / 240.0) * ihy;
    tcz = (-1.0 / 240.0) * ihz;

    for (ix = 2; ix < dimx + 2; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < dimy + 2; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            for (iz = 2; iz < dimz + 2; iz++)
            {
                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = -cc * rptr[ixs + iys + iz];
                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -= 
                    fcx * (rptr[ixms + iys + iz] + rptr[ixps + iys + iz]) +
                    fcy * (rptr[ixs + iyms + iz] + rptr[ixs + iyps + iz]) +
                    fcz * (rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    ecxz * (rptr[ixms + iys + (iz - 1)] + rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] + rptr[ixps + iys + (iz + 1)]) +
                    ecyz * (rptr[ixs + iyms + (iz - 1)] + rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyps + (iz + 1)]) +
                    ecxy * (rptr[ixms + iyms + iz] + rptr[ixms + iyps + iz] +
                            rptr[ixps + iyms + iz] + rptr[ixps + iyps + iz]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    cor * (rptr[ixms + iyms + (iz - 1)] + rptr[ixps + iyms + (iz - 1)] +
                           rptr[ixms + iyps + (iz - 1)] + rptr[ixps + iyps + (iz - 1)] +
                           rptr[ixms + iyms + (iz + 1)] + rptr[ixps + iyms + (iz + 1)] +
                           rptr[ixms + iyps + (iz + 1)] + rptr[ixps + iyps + (iz + 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz]) +
                    fc2y * (rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                    fc2z * (rptr[ixs + iys + (iz - 2)] + rptr[ixs + iys + (iz + 2)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    tcx * (rptr[ixps + iypps + iz] + rptr[ixps + iymms + iz] +
                           rptr[ixms + iypps + iz] + rptr[ixms + iymms + iz] +
                           rptr[ixps + iys + (iz + 2)] + rptr[ixps + iys + (iz - 2)] +
                           rptr[ixms + iys + (iz + 2)] + rptr[ixms + iys + (iz - 2)]) +
                    tcy * (rptr[ixpps + iyps + iz] + rptr[ixmms + iyps + iz] +
                           rptr[ixpps + iyms + iz] + rptr[ixmms + iyms + iz] +
                           rptr[ixs + iyps + (iz + 2)] + rptr[ixs + iyps + (iz - 2)] +
                           rptr[ixs + iyms + (iz + 2)] + rptr[ixs + iyms + (iz - 2)]) +
                    tcz * (rptr[ixpps + iys + (iz + 1)] + rptr[ixmms + iys + (iz + 1)] +
                           rptr[ixpps + iys + (iz - 1)] + rptr[ixmms + iys + (iz - 1)] +
                           rptr[ixs + iypps + (iz + 1)] + rptr[ixs + iymms + (iz + 1)] +
                           rptr[ixs + iypps + (iz - 1)] + rptr[ixs + iymms + (iz - 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c000*rptr[ixs + iys + iz]        *vtot_eig_s[ixs + iys + iz]+
                  2.0*c100*(rptr[ixs + iys + (iz - 1)] *vtot_eig_s[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] *vtot_eig_s[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz]      *vtot_eig_s[ixms + iys + iz] +
                            rptr[ixps + iys + iz]      *vtot_eig_s[ixps + iys + iz] +
                            rptr[ixs + iyms + iz]      *vtot_eig_s[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]      *vtot_eig_s[ixs + iyps + iz]); 

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c110*(rptr[ixps + iyps + iz]       *vtot_eig_s[ixps + iyps + iz] +      
                            rptr[ixps + iyms + iz]       *vtot_eig_s[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz]       *vtot_eig_s[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz]       *vtot_eig_s[ixms + iyms + iz] +
                            rptr[ixps + iys + (iz + 1)]  *vtot_eig_s[ixps + iys + (iz + 1)] +
                            rptr[ixps + iys + (iz - 1)]  *vtot_eig_s[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)]  *vtot_eig_s[ixms + iys + (iz + 1)] +
                            rptr[ixms + iys + (iz - 1)]  *vtot_eig_s[ixms + iys + (iz - 1)] +
                            rptr[ixs + iyps + (iz + 1)]  *vtot_eig_s[ixs + iyps + (iz + 1)] +
                            rptr[ixs + iyps + (iz - 1)]  *vtot_eig_s[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)]  *vtot_eig_s[ixs + iyms + (iz + 1)] + 
                            rptr[ixs + iyms + (iz - 1)]  *vtot_eig_s[ixs + iyms + (iz - 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c200*(rptr[ixs + iys + (iz - 2)] *vtot_eig_s[ixs + iys + (iz - 2)] +
                            rptr[ixs + iys + (iz + 2)] *vtot_eig_s[ixs + iys + (iz + 2)] +
                            rptr[ixmms + iys + iz]     *vtot_eig_s[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz]     *vtot_eig_s[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz]     *vtot_eig_s[ixs + iymms + iz] + 
                            rptr[ixs + iypps + iz]     *vtot_eig_s[ixs + iypps + iz]);


                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c110 * (rptr[ixps + iyps + iz] +
                            rptr[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz] +
                            rptr[ixps + iys + (iz + 1)] +
                            rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] +
                            rptr[ixms + iys + (iz - 1)] +
                            rptr[ixs + iyps + (iz + 1)] +
                            rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyms + (iz - 1)]);

                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (rptr[ixs + iys + (iz - 2)] +
                            rptr[ixs + iys + (iz + 2)] +
                            rptr[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    my_free(rptr);
    return cc;

}

