#include <complex>
#include "const.h"
#include "common_prototypes.h"
#include "FiniteDiff.h"
#include "TradeImages.h"
#include "transition.h"

template double AppCilrSixth<double>(double *, double *, double *, double *, int, int, int, double, double, double);
template double AppCilrSixth<std::complex<double> >(std::complex<double> *, std::complex<double> *, std::complex<double> *, double *, int, int, int, double, double, double);

template  <typename OrbitalType> double AppCilrSixth (OrbitalType *psi, OrbitalType *a_psi, OrbitalType *b_psi, double *vtot, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    double ecxy, ecxz, ecyz, cc, fcx, fcy, fcz, cor;
    double fc2x, fc2y, fc2z, tcx, tcy, tcz;
    double ihx, ihy, ihz;
    double c000, c100, c110, c200;


    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;


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
                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = -cc * psi[ixs + iys + iz];
                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -= 
                    fcx * (psi[ixms + iys + iz] + psi[ixps + iys + iz]) +
                    fcy * (psi[ixs + iyms + iz] + psi[ixs + iyps + iz]) +
                    fcz * (psi[ixs + iys + (iz - 1)] + psi[ixs + iys + (iz + 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    ecxz * (psi[ixms + iys + (iz - 1)] + psi[ixps + iys + (iz - 1)] +
                            psi[ixms + iys + (iz + 1)] + psi[ixps + iys + (iz + 1)]) +
                    ecyz * (psi[ixs + iyms + (iz - 1)] + psi[ixs + iyps + (iz - 1)] +
                            psi[ixs + iyms + (iz + 1)] + psi[ixs + iyps + (iz + 1)]) +
                    ecxy * (psi[ixms + iyms + iz] + psi[ixms + iyps + iz] +
                            psi[ixps + iyms + iz] + psi[ixps + iyps + iz]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    cor * (psi[ixms + iyms + (iz - 1)] + psi[ixps + iyms + (iz - 1)] +
                           psi[ixms + iyps + (iz - 1)] + psi[ixps + iyps + (iz - 1)] +
                           psi[ixms + iyms + (iz + 1)] + psi[ixps + iyms + (iz + 1)] +
                           psi[ixms + iyps + (iz + 1)] + psi[ixps + iyps + (iz + 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    fc2x * (psi[ixmms + iys + iz] + psi[ixpps + iys + iz]) +
                    fc2y * (psi[ixs + iymms + iz] + psi[ixs + iypps + iz]) +
                    fc2z * (psi[ixs + iys + (iz - 2)] + psi[ixs + iys + (iz + 2)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] -=
                    tcx * (psi[ixps + iypps + iz] + psi[ixps + iymms + iz] +
                           psi[ixms + iypps + iz] + psi[ixms + iymms + iz] +
                           psi[ixps + iys + (iz + 2)] + psi[ixps + iys + (iz - 2)] +
                           psi[ixms + iys + (iz + 2)] + psi[ixms + iys + (iz - 2)]) +
                    tcy * (psi[ixpps + iyps + iz] + psi[ixmms + iyps + iz] +
                           psi[ixpps + iyms + iz] + psi[ixmms + iyms + iz] +
                           psi[ixs + iyps + (iz + 2)] + psi[ixs + iyps + (iz - 2)] +
                           psi[ixs + iyms + (iz + 2)] + psi[ixs + iyms + (iz - 2)]) +
                    tcz * (psi[ixpps + iys + (iz + 1)] + psi[ixmms + iys + (iz + 1)] +
                           psi[ixpps + iys + (iz - 1)] + psi[ixmms + iys + (iz - 1)] +
                           psi[ixs + iypps + (iz + 1)] + psi[ixs + iymms + (iz + 1)] +
                           psi[ixs + iypps + (iz - 1)] + psi[ixs + iymms + (iz - 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c000*psi[ixs + iys + iz]        *vtot[ixs + iys + iz]+
                  2.0*c100*(psi[ixs + iys + (iz - 1)] *vtot[ixs + iys + (iz - 1)] +
                            psi[ixs + iys + (iz + 1)] *vtot[ixs + iys + (iz + 1)] +
                            psi[ixms + iys + iz]      *vtot[ixms + iys + iz] +
                            psi[ixps + iys + iz]      *vtot[ixps + iys + iz] +
                            psi[ixs + iyms + iz]      *vtot[ixs + iyms + iz] +
                            psi[ixs + iyps + iz]      *vtot[ixs + iyps + iz]); 

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c110*(psi[ixps + iyps + iz]       *vtot[ixps + iyps + iz] +      
                            psi[ixps + iyms + iz]       *vtot[ixps + iyms + iz] +
                            psi[ixms + iyps + iz]       *vtot[ixms + iyps + iz] +
                            psi[ixms + iyms + iz]       *vtot[ixms + iyms + iz] +
                            psi[ixps + iys + (iz + 1)]  *vtot[ixps + iys + (iz + 1)] +
                            psi[ixps + iys + (iz - 1)]  *vtot[ixps + iys + (iz - 1)] +
                            psi[ixms + iys + (iz + 1)]  *vtot[ixms + iys + (iz + 1)] +
                            psi[ixms + iys + (iz - 1)]  *vtot[ixms + iys + (iz - 1)] +
                            psi[ixs + iyps + (iz + 1)]  *vtot[ixs + iyps + (iz + 1)] +
                            psi[ixs + iyps + (iz - 1)]  *vtot[ixs + iyps + (iz - 1)] +
                            psi[ixs + iyms + (iz + 1)]  *vtot[ixs + iyms + (iz + 1)] + 
                            psi[ixs + iyms + (iz - 1)]  *vtot[ixs + iyms + (iz - 1)]);

                a_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                  2.0*c200*(psi[ixs + iys + (iz - 2)] *vtot[ixs + iys + (iz - 2)] +
                            psi[ixs + iys + (iz + 2)] *vtot[ixs + iys + (iz + 2)] +
                            psi[ixmms + iys + iz]     *vtot[ixmms + iys + iz] +
                            psi[ixpps + iys + iz]     *vtot[ixpps + iys + iz] +
                            psi[ixs + iymms + iz]     *vtot[ixs + iymms + iz] + 
                            psi[ixs + iypps + iz]     *vtot[ixs + iypps + iz]);


                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (psi[ixs + iys + (iz - 1)] +
                            psi[ixs + iys + (iz + 1)] +
                            psi[ixms + iys + iz] +
                            psi[ixps + iys + iz] +
                            psi[ixs + iyms + iz] +
                            psi[ixs + iyps + iz]) + c000 * psi[ixs + iys + iz];

                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c110 * (psi[ixps + iyps + iz] +
                            psi[ixps + iyms + iz] +
                            psi[ixms + iyps + iz] +
                            psi[ixms + iyms + iz] +
                            psi[ixps + iys + (iz + 1)] +
                            psi[ixps + iys + (iz - 1)] +
                            psi[ixms + iys + (iz + 1)] +
                            psi[ixms + iys + (iz - 1)] +
                            psi[ixs + iyps + (iz + 1)] +
                            psi[ixs + iyps + (iz - 1)] +
                            psi[ixs + iyms + (iz + 1)] + psi[ixs + iyms + (iz - 1)]);

                b_psi[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (psi[ixs + iys + (iz - 2)] +
                            psi[ixs + iys + (iz + 2)] +
                            psi[ixmms + iys + iz] +
                            psi[ixpps + iys + iz] +
                            psi[ixs + iymms + iz] + psi[ixs + iypps + iz]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    return cc;

}

