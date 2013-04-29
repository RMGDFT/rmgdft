/************************** SVN Revision Information **************************
 **    $Id: app_cil_fourth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "hybrid.h"


REAL app_cilr_fourth (REAL *psi, REAL *a_psi, REAL *b_psi, REAL *vtot_eig_s, int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz)
{

    int  numgrid, used_alloc=FALSE;
    rmg_double_t *rptr=NULL;

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    REAL ecxy, ecxz, ecyz, cc, fcx, fcy, fcz;
    REAL ihx, ihy, ihz, a1, a2, a3;
    REAL c000, c100;
    int tid;

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case
#else
    tid = 0;
#endif



    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);


    // If rptr is null then we must allocate it here
    if(rptr == NULL) {
        my_malloc (rptr, sbasis + 64, rmg_double_t);
        used_alloc = TRUE;
    }


    if((ct.ibrav != CUBIC_PRIMITIVE) && (ct.ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cil_fourth.\n");
    }


    trade_imagesx (psi, rptr, dimx, dimy, dimz, 1, FULL_FD);


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
    ihy = 1.0 / (gridhy * gridhy * ct.yside * ct.yside);
    ihz = 1.0 / (gridhz * gridhz * ct.zside * ct.zside);

    c000 = 0.5;
    c100 = 1.0 / 12.0;

    if (ct.anisotropy < 1.000001)
    {

        ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
        cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
        fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
        ecxy = (1.0 / 12.0) * (ihx + ihx);
        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;

        for (ix = 1; ix <= dimx; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            for (iy = 1; iy <= dimy; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;

                for (iz = 1; iz <= dimz; iz++)
                {

                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        -cc * rptr[ixs + iys + iz] +
                        -fcx * (rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] +
                                rptr[ixs + iys + (iz - 1)] + 
                                rptr[ixs + iys + (iz + 1)]);

                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        -ecxy * (rptr[ixms + iys + iz - 1] +
                                rptr[ixps + iys + iz - 1] +
                                rptr[ixs + iyms + iz - 1] +
                                rptr[ixs + iyps + iz - 1] +
                                rptr[ixms + iyms + iz] +
                                rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] +
                                rptr[ixps + iyps + iz] +
                                rptr[ixms + iys + iz + 1] +
                                rptr[ixps + iys + iz + 1] +
                                rptr[ixs + iyms + iz + 1] + 
                                rptr[ixs + iyps + iz + 1]);
                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        2.0*c100 * (rptr[ixs + iys + (iz - 1)] *vtot_eig_s[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] *vtot_eig_s[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] *     vtot_eig_s[ixms + iys + iz] +
                                rptr[ixps + iys + iz] *     vtot_eig_s[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] *     vtot_eig_s[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] *     vtot_eig_s[ixs + iyps + iz])+ 
                    2.0*c000 *  rptr[ixs + iys + iz]   *    vtot_eig_s[ixs + iys + iz];

                    b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        c100 * (rptr[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz]) + 
                        c000 *  rptr[ixs + iys + iz];


                }           /* end for */

            }               /* end for */

        }                   /* end for */
    }
    else
    {

        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
        ihy = 1.0 / (gridhy * gridhy * ct.yside * ct.yside);
        ihz = 1.0 / (gridhz * gridhz * ct.zside * ct.zside);

        cc = (-4.0 / 3.0) * (ihx + ihy + ihz);

        fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
        fcy = (5.0 / 6.0) * ihy + (cc / 8.0);
        fcz = (5.0 / 6.0) * ihz + (cc / 8.0);

        ecxy = (1.0 / 12.0) * (ihx + ihy);
        ecxz = (1.0 / 12.0) * (ihx + ihz);
        ecyz = (1.0 / 12.0) * (ihy + ihz);


        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;



        for (ix = 1; ix <= dimx; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            for (iy = 1; iy <= dimy; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;

                for (iz = 1; iz <= dimz; iz++)
                {

                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        -cc * rptr[ixs + iys + iz] ;
                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                        fcx * rptr[ixms + iys + iz] +
                        fcx * rptr[ixps + iys + iz] +
                        fcy * rptr[ixs + iyms + iz] +
                        fcy * rptr[ixs + iyps + iz] +
                        fcz * rptr[ixs + iys + (iz - 1)] + 
                        fcz * rptr[ixs + iys + (iz + 1)];

                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                        ecxz * rptr[ixms + iys + iz - 1] +
                        ecxz * rptr[ixps + iys + iz - 1] +
                        ecyz * rptr[ixs + iyms + iz - 1] +
                        ecyz * rptr[ixs + iyps + iz - 1] +
                        ecxy * rptr[ixms + iyms + iz] +
                        ecxy * rptr[ixms + iyps + iz] +
                        ecxy * rptr[ixps + iyms + iz] +
                        ecxy * rptr[ixps + iyps + iz] +
                        ecxz * rptr[ixms + iys + iz + 1] +
                        ecxz * rptr[ixps + iys + iz + 1] +
                        ecyz * rptr[ixs + iyms + iz + 1] + 
                        ecyz * rptr[ixs + iyps + iz + 1];
                    a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        2.0*c100 * (rptr[ixs + iys + (iz - 1)] *vtot_eig_s[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] *vtot_eig_s[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] *     vtot_eig_s[ixms + iys + iz] +
                                rptr[ixps + iys + iz] *     vtot_eig_s[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] *     vtot_eig_s[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] *     vtot_eig_s[ixs + iyps + iz]) +
                    2.0*c000 *  rptr[ixs + iys + iz]   *    vtot_eig_s[ixs + iys + iz];

                    b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        c100 * (rptr[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz]) + 
                        c000 *  rptr[ixs + iys + iz];


                }           /* end for */

            }               /* end for */

        }                   /* end for */

    }                       /* end if */


    if(used_alloc)
        my_free(rptr);
    return cc;

}                               /* end app_cil */

