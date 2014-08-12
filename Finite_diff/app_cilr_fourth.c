/************************** SVN Revision Information **************************
 **    $Id: app_cil_fourth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "hybrid.h"
#include "FiniteDiff.h"
#include "TradeImages.h"


double app_cilr_fourth (double *psi, double *a_psi, double *b_psi, double *vtot_eig_s, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    int  used_alloc=FALSE, ibrav;
    double *rptr=NULL;

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    double ecxy, ecxz, ecyz, cc, fcx, fcy, fcz;
    double ihx, ihy, ihz, a1, a2, a3, Bz;
    double c000 = 0.5;
    double c100 = 1.0 / 12.0;
    double Bc = 2.0 / 3.0;
    double Bf = 1.0 / 36.0;
    int tid;

    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case

    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);


    ibrav = get_ibrav_type();

    // If rptr is null then we must allocate it here
    if(rptr == NULL) {
        my_malloc (rptr, sbasis + 64, double);
        used_alloc = TRUE;
    }

    trade_imagesx (psi, rptr, dimx, dimy, dimz, 1, FULL_TRADE);


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());


    switch(ibrav) {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
            if (get_anisotropy() < 1.000001)
            {

                cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
                fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
                ecxy = (1.0 / 12.0) * (ihx + ihx);

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
                cc = (-4.0 / 3.0) * (ihx + ihy + ihz);

                fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
                fcy = (5.0 / 6.0) * ihy + (cc / 8.0);
                fcz = (5.0 / 6.0) * ihz + (cc / 8.0);

                ecxy = (1.0 / 12.0) * (ihx + ihy);
                ecxz = (1.0 / 12.0) * (ihx + ihz);
                ecyz = (1.0 / 12.0) * (ihy + ihz);

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
                                2.0*c100 * (
                                        rptr[ixs + iys + (iz - 1)] *vtot_eig_s[ixs + iys + (iz - 1)] +
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
            break;

        case HEXAGONAL:

            Bc = 7.0 / 12.0;
            Bf = 1.0 / 24.0;
            Bz = 1.0 / 12.0;
            cc = ((-3.0 / 4.0) * ihz) - ((5.0 / 3.0) * ihx);
            a1 = ((3.0 / 8.0) * ihz) - ((1.0 / 6.0) * ihx);
            a2 = ((5.0 / 18.0) * ihx) - ((1.0 / 24.0) * ihz);
            a3 = ((1.0 / 48.0) * ihz) + ((1.0 / 36.0) * ihx);
            cc *= 2.0;
            a1 *= 2.0;
            a2 *= 2.0;
            a3 *= 2.0;

            for(ix = 1;ix <= dimx;ix++) {
              ixs = ix * incx;
              ixms = (ix - 1) * incx;
              ixps = (ix + 1) * incx;
              for(iy = 1;iy <= dimy;iy++) {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                for(iz = 1;iz <= dimz;iz++) {
  
                    b_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] = 
          
                             Bc * rptr[ixs + iys + iz] +
                             Bz * rptr[ixs + iys + iz - 1] +
                             Bz * rptr[ixs + iys + iz + 1] +
                             Bf * rptr[ixps + iys + iz] +
                             Bf * rptr[ixps + iyms + iz] +
                             Bf * rptr[ixs + iyms + iz] +
                             Bf * rptr[ixms + iys + iz] +
                             Bf * rptr[ixms + iyps + iz] +
                             Bf * rptr[ixs + iyps + iz];


                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] =
  
                      -cc * rptr[ixs + iys + iz];

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=

                         a3 * (rptr[ixps + iys + iz-1] +
                               rptr[ixps + iyms + iz-1] +
                               rptr[ixs + iyms + iz-1] +
                               rptr[ixms + iys + iz-1] +
                               rptr[ixms + iyps + iz-1] +
                               rptr[ixs + iyps + iz-1] +

                               rptr[ixps + iys + iz+1] +
                               rptr[ixps + iyms + iz+1] +
                               rptr[ixs + iyms + iz+1] +
                               rptr[ixms + iys + iz+1] +
                               rptr[ixms + iyps + iz+1] +
                               rptr[ixs + iyps + iz+1]);
     
                                                            
                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=
                         a2 * (rptr[ixps + iys + iz] +
                               rptr[ixps + iyms + iz] +
                               rptr[ixs + iyms + iz] +
                               rptr[ixms + iys + iz] +
                               rptr[ixms + iyps + iz] +
                               rptr[ixs + iyps + iz]);

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=
                         a1 * (rptr[ixs + iys + iz-1] + rptr[ixs + iys + iz+1]);

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] +=
                             2.0 * Bc * rptr[ixs + iys + iz] * vtot_eig_s[ixs + iys + iz] +
                             2.0 * Bz * rptr[ixs + iys + iz - 1] * vtot_eig_s[ixs + iys + iz - 1] +
                             2.0 * Bz * rptr[ixs + iys + iz + 1] * vtot_eig_s[ixs + iys + iz + 1] +
                             2.0 * Bf * rptr[ixps + iys + iz] * vtot_eig_s[ixps + iys + iz] +
                             2.0 * Bf * rptr[ixps + iyms + iz] * vtot_eig_s[ixps + iyms + iz] +
                             2.0 * Bf * rptr[ixs + iyms + iz] * vtot_eig_s[ixs + iyms + iz] +
                             2.0 * Bf * rptr[ixms + iys + iz] * vtot_eig_s[ixms + iys + iz] +
                             2.0 * Bf * rptr[ixms + iyps + iz] * vtot_eig_s[ixms + iyps + iz] +
                             2.0 * Bf * rptr[ixs + iyps + iz] * vtot_eig_s[ixs + iyps + iz];
                                                            
                } /* end for */
  
              } /* end for */

            } /* end for */
            break;

        case CUBIC_FC:

            /* Compute coefficients for this grid spacing */
            cc = (-34.0 / 6.0) * ihx;
            a1 = (4.0 / 9.0) * ihx;
            a2 = (1.0 / 18.0) * ihx;


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
                            -cc * rptr[ixs + iys + iz];

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                            a1 * (rptr[ixms + iys + iz] +
                                  rptr[ixms + iys + iz + 1] +
                                  rptr[ixms + iyps + iz] +
                                  rptr[ixs + iyms + iz] +
                                  rptr[ixs + iyms + iz + 1] +
                                  rptr[ixs + iys + iz - 1] +
                                  rptr[ixs + iys + iz + 1] +
                                  rptr[ixs + iyps + iz - 1] +
                                  rptr[ixs + iyps + iz] +
                                  rptr[ixps + iyms + iz] +
                                  rptr[ixps + iy * incy + iz - 1] + 
                                  rptr[ixps + iy * incy + iz]);

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                            a2 * (rptr[ixms + iyms + iz + 1] +
                                  rptr[ixms + iyps + iz - 1] +
                                  rptr[ixms + iyps + iz + 1] +
                                  rptr[ixps + iyms + iz - 1] +
                                  rptr[ixps + iyms + iz + 1] + 
                                  rptr[ixps + iyps + iz - 1]);

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                2.0 * Bc * rptr[ixs + iys + iz] * vtot_eig_s[ixs + iys + iz] +
                                2.0 * Bf * 
                                     (rptr[ixms + iys + iz] * vtot_eig_s[ixms + iys + iz] +
                                      rptr[ixms + iys + iz + 1] * vtot_eig_s[ixms + iys + iz + 1] +
                                      rptr[ixms + iyps + iz] * vtot_eig_s[ixms + iyps + iz] +
                                      rptr[ixs + iyms + iz] * vtot_eig_s[ixs + iyms + iz] +
                                      rptr[ixs + iyms + iz + 1] * vtot_eig_s[ixs + iyms + iz + 1] +
                                      rptr[ixs + iys + iz - 1] * vtot_eig_s[ixs + iys + iz - 1] +
                                      rptr[ixs + iys + iz + 1] * vtot_eig_s[ixs + iys + iz + 1] +
                                      rptr[ixs + iyps + iz - 1] * vtot_eig_s[ixs + iyps + iz - 1] +
                                      rptr[ixs + iyps + iz] * vtot_eig_s[ixs + iyps + iz] +
                                      rptr[ixps + iyms + iz] * vtot_eig_s[ixps + iyms + iz] +
                                      rptr[ixps + iys + iz - 1] * vtot_eig_s[ixps + iys + iz - 1] +
                                      rptr[ixps + iys + iz] * vtot_eig_s[ixps + iys + iz]);

                        b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                Bc * rptr[ixs + iys + iz] +
                                Bf * (rptr[ixms + iys + iz] +
                                      rptr[ixms + iys + iz + 1] +
                                      rptr[ixms + iyps + iz] +
                                      rptr[ixs + iyms + iz] +
                                      rptr[ixs + iyms + iz + 1] +
                                      rptr[ixs + iys + iz - 1] +
                                      rptr[ixs + iys + iz + 1] +
                                      rptr[ixs + iyps + iz - 1] +
                                      rptr[ixs + iyps + iz] +
                                      rptr[ixps + iyms + iz] +
                                      rptr[ixps + iys + iz - 1] +
                                      rptr[ixps + iys + iz]);


                    }           /* end for */

                }               /* end for */

            }                   /* end for */

            break;

        default:
            error_handler("Grid symmetry not programmed yet in app_cil_fourth.\n");

    } // end switch
    if(used_alloc)
        my_free(rptr);
    return cc;

}                               /* end app_cil */

