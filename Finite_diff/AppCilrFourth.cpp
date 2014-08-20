#include <complex>
#include "const.h"
#include "common_prototypes.h"
#include "FiniteDiff.h"
#include "TradeImages.h"
#include "transition.h"

template double AppCilrFourth<double>(double *, double *, double *, double *, int, int, int, double, double, double);
template double AppCilrFourth<std::complex<double> >(std::complex<double> *, std::complex<double> *, std::complex<double> *, double *, int, int, int, double, double, double);

template  <typename OrbitalType> double AppCilrFourth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    double cc = 0.0;
    double ecxy, ecxz, ecyz, fcx, fcy, fcz;
    double ihx, ihy, ihz, a1, a2, a3, Bz;
    double c000 = 0.5;
    double c100 = 1.0 / 12.0;
    double Bc = 2.0 / 3.0;
    double Bf = 1.0 / 36.0;

    int ibrav = get_ibrav_type();


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
                                -cc * psi[ixs + iys + iz] +
                                -fcx * (psi[ixms + iys + iz] +
                                        psi[ixps + iys + iz] +
                                        psi[ixs + iyms + iz] +
                                        psi[ixs + iyps + iz] +
                                        psi[ixs + iys + (iz - 1)] + 
                                        psi[ixs + iys + (iz + 1)]);

                            a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                -ecxy * (psi[ixms + iys + iz - 1] +
                                        psi[ixps + iys + iz - 1] +
                                        psi[ixs + iyms + iz - 1] +
                                        psi[ixs + iyps + iz - 1] +
                                        psi[ixms + iyms + iz] +
                                        psi[ixms + iyps + iz] +
                                        psi[ixps + iyms + iz] +
                                        psi[ixps + iyps + iz] +
                                        psi[ixms + iys + iz + 1] +
                                        psi[ixps + iys + iz + 1] +
                                        psi[ixs + iyms + iz + 1] + 
                                        psi[ixs + iyps + iz + 1]);
                            a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                2.0*c100 * (psi[ixs + iys + (iz - 1)] *vtot[ixs + iys + (iz - 1)] +
                                        psi[ixs + iys + (iz + 1)] *vtot[ixs + iys + (iz + 1)] +
                                        psi[ixms + iys + iz] *     vtot[ixms + iys + iz] +
                                        psi[ixps + iys + iz] *     vtot[ixps + iys + iz] +
                                        psi[ixs + iyms + iz] *     vtot[ixs + iyms + iz] +
                                        psi[ixs + iyps + iz] *     vtot[ixs + iyps + iz])+ 
                            2.0*c000 *  psi[ixs + iys + iz]   *    vtot[ixs + iys + iz];

                            b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                c100 * (psi[ixs + iys + (iz - 1)] +
                                        psi[ixs + iys + (iz + 1)] +
                                        psi[ixms + iys + iz] +
                                        psi[ixps + iys + iz] +
                                        psi[ixs + iyms + iz] +
                                        psi[ixs + iyps + iz]) + 
                                c000 *  psi[ixs + iys + iz];


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
                                -cc * psi[ixs + iys + iz] ;
                            a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                                fcx * psi[ixms + iys + iz] +
                                fcx * psi[ixps + iys + iz] +
                                fcy * psi[ixs + iyms + iz] +
                                fcy * psi[ixs + iyps + iz] +
                                fcz * psi[ixs + iys + (iz - 1)] + 
                                fcz * psi[ixs + iys + (iz + 1)];

                            a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                                ecxz * psi[ixms + iys + iz - 1] +
                                ecxz * psi[ixps + iys + iz - 1] +
                                ecyz * psi[ixs + iyms + iz - 1] +
                                ecyz * psi[ixs + iyps + iz - 1] +
                                ecxy * psi[ixms + iyms + iz] +
                                ecxy * psi[ixms + iyps + iz] +
                                ecxy * psi[ixps + iyms + iz] +
                                ecxy * psi[ixps + iyps + iz] +
                                ecxz * psi[ixms + iys + iz + 1] +
                                ecxz * psi[ixps + iys + iz + 1] +
                                ecyz * psi[ixs + iyms + iz + 1] + 
                                ecyz * psi[ixs + iyps + iz + 1];

                            a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                2.0*c100 * (
                                        psi[ixs + iys + (iz - 1)] *vtot[ixs + iys + (iz - 1)] +
                                        psi[ixs + iys + (iz + 1)] *vtot[ixs + iys + (iz + 1)] +
                                        psi[ixms + iys + iz] *     vtot[ixms + iys + iz] +
                                        psi[ixps + iys + iz] *     vtot[ixps + iys + iz] +
                                        psi[ixs + iyms + iz] *     vtot[ixs + iyms + iz] +
                                        psi[ixs + iyps + iz] *     vtot[ixs + iyps + iz]) +
                            2.0*c000 *  psi[ixs + iys + iz]   *    vtot[ixs + iys + iz];

                            b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                c100 * (psi[ixs + iys + (iz - 1)] +
                                        psi[ixs + iys + (iz + 1)] +
                                        psi[ixms + iys + iz] +
                                        psi[ixps + iys + iz] +
                                        psi[ixs + iyms + iz] +
                                        psi[ixs + iyps + iz]) + 
                                c000 *  psi[ixs + iys + iz];


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
          
                             Bc * psi[ixs + iys + iz] +
                             Bz * psi[ixs + iys + iz - 1] +
                             Bz * psi[ixs + iys + iz + 1] +
                             Bf * psi[ixps + iys + iz] +
                             Bf * psi[ixps + iyms + iz] +
                             Bf * psi[ixs + iyms + iz] +
                             Bf * psi[ixms + iys + iz] +
                             Bf * psi[ixms + iyps + iz] +
                             Bf * psi[ixs + iyps + iz];


                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] =
  
                      -cc * psi[ixs + iys + iz];

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=

                         a3 * (psi[ixps + iys + iz-1] +
                               psi[ixps + iyms + iz-1] +
                               psi[ixs + iyms + iz-1] +
                               psi[ixms + iys + iz-1] +
                               psi[ixms + iyps + iz-1] +
                               psi[ixs + iyps + iz-1] +

                               psi[ixps + iys + iz+1] +
                               psi[ixps + iyms + iz+1] +
                               psi[ixs + iyms + iz+1] +
                               psi[ixms + iys + iz+1] +
                               psi[ixms + iyps + iz+1] +
                               psi[ixs + iyps + iz+1]);
     
                                                            
                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=
                         a2 * (psi[ixps + iys + iz] +
                               psi[ixps + iyms + iz] +
                               psi[ixs + iyms + iz] +
                               psi[ixms + iys + iz] +
                               psi[ixms + iyps + iz] +
                               psi[ixs + iyps + iz]);

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] -=
                         a1 * (psi[ixs + iys + iz-1] + psi[ixs + iys + iz+1]);

                    a_psi[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] +=
                             2.0 * Bc * psi[ixs + iys + iz] * vtot[ixs + iys + iz] +
                             2.0 * Bz * psi[ixs + iys + iz - 1] * vtot[ixs + iys + iz - 1] +
                             2.0 * Bz * psi[ixs + iys + iz + 1] * vtot[ixs + iys + iz + 1] +
                             2.0 * Bf * psi[ixps + iys + iz] * vtot[ixps + iys + iz] +
                             2.0 * Bf * psi[ixps + iyms + iz] * vtot[ixps + iyms + iz] +
                             2.0 * Bf * psi[ixs + iyms + iz] * vtot[ixs + iyms + iz] +
                             2.0 * Bf * psi[ixms + iys + iz] * vtot[ixms + iys + iz] +
                             2.0 * Bf * psi[ixms + iyps + iz] * vtot[ixms + iyps + iz] +
                             2.0 * Bf * psi[ixs + iyps + iz] * vtot[ixs + iyps + iz];
                                                            
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
                            -cc * psi[ixs + iys + iz];

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                            a1 * (psi[ixms + iys + iz] +
                                  psi[ixms + iys + iz + 1] +
                                  psi[ixms + iyps + iz] +
                                  psi[ixs + iyms + iz] +
                                  psi[ixs + iyms + iz + 1] +
                                  psi[ixs + iys + iz - 1] +
                                  psi[ixs + iys + iz + 1] +
                                  psi[ixs + iyps + iz - 1] +
                                  psi[ixs + iyps + iz] +
                                  psi[ixps + iyms + iz] +
                                  psi[ixps + iy * incy + iz - 1] + 
                                  psi[ixps + iy * incy + iz]);

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] -=
                            a2 * (psi[ixms + iyms + iz + 1] +
                                  psi[ixms + iyps + iz - 1] +
                                  psi[ixms + iyps + iz + 1] +
                                  psi[ixps + iyms + iz - 1] +
                                  psi[ixps + iyms + iz + 1] + 
                                  psi[ixps + iyps + iz - 1]);

                        a_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                2.0 * Bc * psi[ixs + iys + iz] * vtot[ixs + iys + iz] +
                                2.0 * Bf * 
                                     (psi[ixms + iys + iz] * vtot[ixms + iys + iz] +
                                      psi[ixms + iys + iz + 1] * vtot[ixms + iys + iz + 1] +
                                      psi[ixms + iyps + iz] * vtot[ixms + iyps + iz] +
                                      psi[ixs + iyms + iz] * vtot[ixs + iyms + iz] +
                                      psi[ixs + iyms + iz + 1] * vtot[ixs + iyms + iz + 1] +
                                      psi[ixs + iys + iz - 1] * vtot[ixs + iys + iz - 1] +
                                      psi[ixs + iys + iz + 1] * vtot[ixs + iys + iz + 1] +
                                      psi[ixs + iyps + iz - 1] * vtot[ixs + iyps + iz - 1] +
                                      psi[ixs + iyps + iz] * vtot[ixs + iyps + iz] +
                                      psi[ixps + iyms + iz] * vtot[ixps + iyms + iz] +
                                      psi[ixps + iys + iz - 1] * vtot[ixps + iys + iz - 1] +
                                      psi[ixps + iys + iz] * vtot[ixps + iys + iz]);

                        b_psi[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                Bc * psi[ixs + iys + iz] +
                                Bf * (psi[ixms + iys + iz] +
                                      psi[ixms + iys + iz + 1] +
                                      psi[ixms + iyps + iz] +
                                      psi[ixs + iyms + iz] +
                                      psi[ixs + iyms + iz + 1] +
                                      psi[ixs + iys + iz - 1] +
                                      psi[ixs + iys + iz + 1] +
                                      psi[ixs + iyps + iz - 1] +
                                      psi[ixs + iyps + iz] +
                                      psi[ixps + iyms + iz] +
                                      psi[ixps + iys + iz - 1] +
                                      psi[ixps + iys + iz]);


                    }           /* end for */

                }               /* end for */

            }                   /* end for */

            break;

        default:
            rmg_error_handler(__FILE__, __LINE__, "Grid symmetry not programmed yet in AppCilrFourth.\n");

    } // end switch
    return cc;

}

