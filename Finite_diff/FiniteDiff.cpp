
#include "fixed_dims.h"
#include "BaseGrid.h"
#include "FiniteDiff.h"
#include <cmath>
#include <complex>

using namespace std;

extern "C" {
  rmg_double_t get_xside(void);
  rmg_double_t get_yside(void);
  rmg_double_t get_zside(void);
  int get_PX0_GRID(void);
  int get_PY0_GRID(void);
  int get_PZ0_GRID(void);
  void rmg_error_handler(const char *message);
}

template <typename RmgType>
rmg_double_t FD_app_cil_sixth_standard (RmgType *rptr, RmgType *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    rmg_double_t ecxy, ecxz, ecyz, cc, fcx, fcy, fcz, cor;
    rmg_double_t fc2x, fc2y, fc2z, tcx, tcy, tcz;
    rmg_double_t ihx, ihy, ihz;

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
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                    fcx * (rptr[ixms + iys + iz] + rptr[ixps + iys + iz]) +
                    fcy * (rptr[ixs + iyms + iz] + rptr[ixs + iyps + iz]) +
                    fcz * (rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    ecxz * (rptr[ixms + iys + (iz - 1)] + rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] + rptr[ixps + iys + (iz + 1)]) +
                    ecyz * (rptr[ixs + iyms + (iz - 1)] + rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyps + (iz + 1)]) +
                    ecxy * (rptr[ixms + iyms + iz] + rptr[ixms + iyps + iz] +
                            rptr[ixps + iyms + iz] + rptr[ixps + iyps + iz]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    cor * (rptr[ixms + iyms + (iz - 1)] + rptr[ixps + iyms + (iz - 1)] +
                           rptr[ixms + iyps + (iz - 1)] + rptr[ixps + iyps + (iz - 1)] +
                           rptr[ixms + iyms + (iz + 1)] + rptr[ixps + iyms + (iz + 1)] +
                           rptr[ixms + iyps + (iz + 1)] + rptr[ixps + iyps + (iz + 1)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz]) +
                    fc2y * (rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                    fc2z * (rptr[ixs + iys + (iz - 2)] + rptr[ixs + iys + (iz + 2)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
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

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    return cc;
}



// Version with loop dimensions set at compile time
template <typename RmgType>
rmg_double_t FD_app_cil_sixth_global (RmgType * rptr, RmgType * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{


    int iz, ix, iy, incx, incy, incxr, incyr, ibrav;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    rmg_double_t ecxy, ecxz, ecyz, cc, fcx, fcy, fcz, cor;
    rmg_double_t fc2x, fc2y, fc2z, tcx, tcy, tcz;
    rmg_double_t ihx, ihy, ihz;
    rmg_double_t rz, rzms, rzps, rzpps;
    rmg_double_t rfc1, rbc1, rbc2, rd1, rd2, rd3, rd4;
    rmg_double_t td1, td2, td3, td4, td5, td6, td7, td8, tdx;

    BaseGrid G;

    ibrav = G.get_ibrav_type();

    incx = (FIXED_ZDIM + 4) * (FIXED_YDIM + 4);
    incy = FIXED_ZDIM + 4;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

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

    // Handle the general case first
    if((FIXED_ZDIM % 4) || (ibrav != CUBIC_PRIMITIVE)) {

        for (ix = 2; ix < FIXED_XDIM + 2; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            ixmms = (ix - 2) * incx;
            ixpps = (ix + 2) * incx;

            for (iy = 2; iy < FIXED_YDIM + 2; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                iymms = (iy - 2) * incy;
                iypps = (iy + 2) * incy;

                for (iz = 2; iz < FIXED_ZDIM + 2; iz++)
                {

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                        fcx * (rptr[ixms + iys + iz] + rptr[ixps + iys + iz]) +
                        fcy * (rptr[ixs + iyms + iz] + rptr[ixs + iyps + iz]) +
                        fcz * (rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        ecxz * (rptr[ixms + iys + (iz - 1)] + rptr[ixps + iys + (iz - 1)] +
                                rptr[ixms + iys + (iz + 1)] + rptr[ixps + iys + (iz + 1)]) +
                        ecyz * (rptr[ixs + iyms + (iz - 1)] + rptr[ixs + iyps + (iz - 1)] +
                                rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyps + (iz + 1)]) +
                        ecxy * (rptr[ixms + iyms + iz] + rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] + rptr[ixps + iyps + iz]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        cor * (rptr[ixms + iyms + (iz - 1)] + rptr[ixps + iyms + (iz - 1)] +
                               rptr[ixms + iyps + (iz - 1)] + rptr[ixps + iyps + (iz - 1)] +
                               rptr[ixms + iyms + (iz + 1)] + rptr[ixps + iyms + (iz + 1)] +
                               rptr[ixms + iyps + (iz + 1)] + rptr[ixps + iyps + (iz + 1)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz]) +
                        fc2y * (rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                        fc2z * (rptr[ixs + iys + (iz - 2)] + rptr[ixs + iys + (iz + 2)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
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

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

        return cc;

    }

    // Optimized case for dimz divisible by 4 and cubic primitive grid

    for (ix = 2; ix < FIXED_XDIM + 2; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < FIXED_YDIM + 2; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            // Compute the middle set of edges (2nd nn) before entry to the loop
            rz = rptr[ixms + iys + 2] +
                 rptr[ixps + iys + 2] +
                 rptr[ixs + iyms + 2] +
                 rptr[ixs + iyps + 2];

            // Compute the trailing set of edges (2nd nn) before entry to the  loop
            rzms = rptr[ixps + iys + 1] +
                   rptr[ixms + iys + 1] +
                   rptr[ixs + iyps + 1] +
                   rptr[ixs + iyms + 1];

            // Compute the 4 nn with same z value as loop index on entry
            rfc1 = rptr[ixms + iyms + 2] + rptr[ixms + iyps + 2] +
                   rptr[ixps + iyms + 2] + rptr[ixps + iyps + 2];

            // Compute a pair trailing sets of corners on entry
            rbc1 = rptr[ixms + iyms + 1] + rptr[ixps + iyms + 1] +
                   rptr[ixms + iyps + 1] + rptr[ixps + iyps + 1];
            rbc2 = rptr[ixms + iyms + 2] + rptr[ixps + iyms + 2] +
                   rptr[ixms + iyps + 2] + rptr[ixps + iyps + 2];

            rd1 = rptr[ixpps + iys + 1] + rptr[ixmms + iys + 1] +
                  rptr[ixs + iypps + 1] + rptr[ixs + iymms + 1];

            rd4 = rptr[ixpps + iys + 2] + rptr[ixmms + iys + 2] +
                  rptr[ixs + iypps + 2] + rptr[ixs + iymms + 2];


            td2 =   rptr[ixps + iys] +
                           rptr[ixms + iys] +
                           rptr[ixs + iyps] +
                           rptr[ixs + iyms];

            td4 =   rptr[ixps + iys + 1] +
                           rptr[ixms + iys + 1] +
                           rptr[ixs + iyps + 1] +
                           rptr[ixs + iyms + 1];

            td6 =   rptr[ixps + iys + 2] +
                           rptr[ixms + iys + 2] +
                           rptr[ixs + iyps + 2] +
                           rptr[ixs + iyms + 2];

            td8 =   rptr[ixps + iys + 3] +
                           rptr[ixms + iys + 3] +
                           rptr[ixs + iyps + 3] +
                           rptr[ixs + iyms + 3];

            td1 =   rptr[ixps + iys + 4] +
                           rptr[ixms + iys + 4] +
                           rptr[ixs + iyps + 4] +
                           rptr[ixs + iyms + 4];

            td3 =   rptr[ixps + iys + 5] +
                           rptr[ixms + iys + 5] +
                           rptr[ixs + iyps + 5] +
                           rptr[ixs + iyms + 5];


            td5 =   rptr[ixps + iys + 6] +
                           rptr[ixms + iys + 6] +
                           rptr[ixs + iyps + 6] +
                           rptr[ixs + iyms + 6];


            td7 =   rptr[ixps + iys + 7] +
                           rptr[ixms + iys + 7] +
                           rptr[ixs + iyps + 7] +
                           rptr[ixs + iyms + 7];



            for (iz = 2; iz < FIXED_ZDIM + 2; iz+=4)
            {

                tdx =   rptr[ixps + iys + iz + 6] +
                        rptr[ixms + iys + iz + 6] +
                        rptr[ixs + iyps + iz + 6] +
                        rptr[ixs + iyms + iz + 6];

                // Forward set of edges
                rzps = rptr[ixps + iys + iz + 1] +
                       rptr[ixms + iys + iz + 1] +
                       rptr[ixs + iyps + iz + 1] +
                       rptr[ixs + iyms + iz + 1];

                // First
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                    fcx * rz +
                    fcx * (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 1]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    ecxz * rzms + ecxz * rzps + ecxz * rfc1;

                // Compute the forward set of corners
                rfc1 = rptr[ixms + iyms + iz + 1] + rptr[ixps + iyms + iz + 1] +
                       rptr[ixms + iyps + iz + 1] + rptr[ixps + iyps + iz + 1];

                rd3 = rptr[ixpps + iys + iz + 1] + rptr[ixmms + iys + iz + 1] +
                          rptr[ixs + iypps + iz + 1] + rptr[ixs + iymms + iz + 1];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    cor * rbc1 + cor * rfc1;
                rbc1 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                    fc2x * (rptr[ixs + iys + iz - 2] + rptr[ixs + iys + iz + 2]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    tcx * rptr[ixps + iypps + iz] + tcx * rptr[ixps + iymms + iz] +
                    tcx * rptr[ixms + iypps + iz] + tcx * rptr[ixms + iymms + iz] +
                    tcx * rptr[ixpps + iyps + iz] + tcx * rptr[ixmms + iyps + iz] +
                    tcx * rptr[ixpps + iyms + iz] + tcx * rptr[ixmms + iyms + iz] +
                    tcx * rd1 +
                    tcx * rd3 +
                    tcx * (td1 + td2);

                td2 = td1;
                td1 = tdx;
                tdx =   rptr[ixps + iys + iz + 7] +
                        rptr[ixms + iys + iz + 7] +
                        rptr[ixs + iyps + iz + 7] +
                        rptr[ixs + iyms + iz + 7];

                // Compute another set of forward edges here
                rzpps = rptr[ixps + iys + iz + 2] +
                        rptr[ixms + iys + iz + 2] +
                        rptr[ixs + iyps + iz + 2] +
                        rptr[ixs + iyms + iz + 2];

                // Second
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] = cc * rptr[ixs + iys + iz + 1] +
                    fcx * rzps +
                    fcx * (rptr[ixs + iys + iz] + rptr[ixs + iys + iz + 2]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    ecxz * rz +  ecxz * rzpps +  ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 2] + rptr[ixps + iyms + iz + 2] +
                       rptr[ixms + iyps + iz + 2] + rptr[ixps + iyps + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    cor * rbc2 + cor * rfc1;
                rbc2 = rfc1;

                rd2 = rptr[ixpps + iys + iz + 2] + rptr[ixmms + iys + iz + 2] +
                      rptr[ixs + iypps + iz + 2] + rptr[ixs + iymms + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    fc2x * (rptr[ixmms + iys + iz + 1] + rptr[ixpps + iys + iz + 1] +
                            rptr[ixs + iymms + iz + 1] + rptr[ixs + iypps + iz + 1]) +
                    fc2x *  (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 3]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    tcx * rptr[ixps + iypps + iz + 1] + tcx * rptr[ixps + iymms + iz + 1] +
                    tcx * rptr[ixms + iypps + iz + 1] + tcx * rptr[ixms + iymms + iz + 1] +
                    tcx * rptr[ixpps + iyps + iz + 1] + tcx * rptr[ixmms + iyps + iz + 1] +
                    tcx * rptr[ixpps + iyms + iz + 1] + tcx * rptr[ixmms + iyms + iz + 1] +
                    tcx * rd4 +
                    tcx * rd2 +
                    tcx * (td3 + td4);

                td4 = td3;
                td3 = tdx;
                tdx =   rptr[ixps + iys + iz + 8] +
                           rptr[ixms + iys + iz + 8] +
                           rptr[ixs + iyps + iz + 8] +
                           rptr[ixs + iyms + iz + 8];

                // Compute the forward set of edges here which becomes the trailing set at the top
                rzms = rptr[ixms + iys + iz+3] +
                       rptr[ixps + iys + iz+3] +
                       rptr[ixs + iyms + iz+3] +
                       rptr[ixs + iyps + iz+3];

                // Third
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] = cc * rptr[ixs + iys + iz + 2] +
                    fcx * rzpps +
                    fcx * (rptr[ixs + iys + iz + 1] + rptr[ixs + iys + iz + 3]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    ecxz * rzps + ecxz * rzms + ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 3] + rptr[ixps + iyms + iz + 3] +
                       rptr[ixms + iyps + iz + 3] + rptr[ixps + iyps + iz + 3];

                rd1 = rptr[ixpps + iys + iz + 3] + rptr[ixmms + iys + iz + 3] +
                      rptr[ixs + iypps + iz + 3] + rptr[ixs + iymms + iz + 3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    cor * rbc1 + cor * rfc1;
                rbc1 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    fc2x * (rptr[ixmms + iys + iz + 2] + rptr[ixpps + iys + iz + 2] +
                            rptr[ixs + iymms + iz + 2] + rptr[ixs + iypps + iz + 2]) +
                    fc2x *  (rptr[ixs + iys + iz] + rptr[ixs + iys + iz + 4]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    tcx * rptr[ixps + iypps + iz + 2] + tcx * rptr[ixps + iymms + iz + 2] +
                    tcx * rptr[ixms + iypps + iz + 2] + tcx * rptr[ixms + iymms + iz + 2] +
                    tcx * rptr[ixpps + iyps + iz + 2] + tcx * rptr[ixmms + iyps + iz + 2] +
                    tcx * rptr[ixpps + iyms + iz + 2] + tcx * rptr[ixmms + iyms + iz + 2] +
                    tcx * rd3 +
                    tcx * rd1 +
                    tcx * (td5 + td6);

                td6 = td5;
                td5 = tdx;
                tdx =   rptr[ixps + iys + iz + 9] +
                        rptr[ixms + iys + iz + 9] +
                        rptr[ixs + iyps + iz + 9] +
                        rptr[ixs + iyms + iz + 9];

                // Compute the forward set of edges here which becomes the middle set at the top
                rz = rptr[ixps + iys + iz + 4] +
                     rptr[ixms + iys + iz + 4] +
                     rptr[ixs + iyps + iz + 4] +
                     rptr[ixs + iyms + iz + 4];

                // Fourth
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] = cc * rptr[ixs + iys + iz + 3] +
                    fcx * rzms +
                    fcx * (rptr[ixs + iys + iz + 2] + rptr[ixs + iys + iz + 4]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    ecxz * rzpps + ecxz * rz + ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 4] + rptr[ixps + iyms + iz + 4] +
                       rptr[ixms + iyps + iz + 4] + rptr[ixps + iyps + iz + 4];

                rd4 = rptr[ixpps + iys + iz + 4] + rptr[ixmms + iys + iz + 4] +
                      rptr[ixs + iypps + iz + 4] + rptr[ixs + iymms + iz + 4];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    cor * rbc2 + cor * rfc1;
                rbc2 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    fc2x * (rptr[ixmms + iys + iz + 3] + rptr[ixpps + iys + iz + 3] +
                            rptr[ixs + iymms + iz + 3] + rptr[ixs + iypps + iz + 3]) +
                    fc2x * (rptr[ixs + iys + iz + 1] + rptr[ixs + iys + iz + 5]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    tcx * rptr[ixps + iypps + iz + 3] + tcx * rptr[ixps + iymms + iz + 3] +
                    tcx * rptr[ixms + iypps + iz + 3] + tcx * rptr[ixms + iymms + iz + 3] +
                    tcx * rptr[ixpps + iyps + iz + 3] + tcx * rptr[ixmms + iyps + iz + 3] +
                    tcx * rptr[ixpps + iyms + iz + 3] + tcx * rptr[ixmms + iyms + iz + 3] +
                    tcx * rd2 +
                    tcx * rd4 +
                    tcx * (td7 + td8);

                td8 = td7;
                td7 = tdx;



            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    return cc;

}



template <typename RmgType>
void FD_app_cir_sixth_standard (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx, ixmms, ixpps, iymms, iypps;
    int incyr, incxr;
    rmg_double_t c000, c100, c110, c200;

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;
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

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
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

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (rptr[ixs + iys + (iz - 2)] +
                            rptr[ixs + iys + (iz + 2)] +
                            rptr[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}


template <typename RmgType>
void FD_app_cir_sixth_global (RmgType * rptr, RmgType * b)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx, ixmms, ixpps, iymms, iypps;
    int incyr, incxr;
    rmg_double_t rz, rzps, rzms, rzpps;
    rmg_double_t c000, c100, c110, c200;

    incx = (FIXED_ZDIM + 4) * (FIXED_YDIM + 4);
    incy = FIXED_ZDIM + 4;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;


    // Handle the general case first
    if(FIXED_ZDIM % 4) {

        for (ix = 2; ix < FIXED_XDIM + 2; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            ixmms = (ix - 2) * incx;
            ixpps = (ix + 2) * incx;

            for (iy = 2; iy < FIXED_YDIM + 2; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                iymms = (iy - 2) * incy;
                iypps = (iy + 2) * incy;

                for (iz = 2; iz < FIXED_ZDIM + 2; iz++)
                {

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                        c100 * (rptr[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
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

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        c200 * (rptr[ixs + iys + (iz - 2)] +
                                rptr[ixs + iys + (iz + 2)] +
                                rptr[ixmms + iys + iz] +
                                rptr[ixpps + iys + iz] +
                                rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]);

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

        return;
    }


    // Optimized case for dimz divisible by 4
    for (ix = 2; ix < FIXED_XDIM + 2; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < FIXED_YDIM + 2; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            // Compute the middle set of edges before entry to the loop
            rz = rptr[ixms + iys + 2] +
                 rptr[ixps + iys + 2] +
                 rptr[ixs + iyms + 2] +
                 rptr[ixs + iyps + 2];

            // Compute the trailing set of edges before entry to the  loop
            rzms = rptr[ixps + iys + 1] +
                   rptr[ixms + iys + 1] +
                   rptr[ixs + iyps + 1] +
                   rptr[ixs + iyms + 1];


            for (iz = 2; iz < FIXED_ZDIM + 2; iz+=4)
            {

                // Forward set of edges
                rzps = rptr[ixps + iys + iz + 1] +
                       rptr[ixms + iys + iz + 1] +
                       rptr[ixs + iyps + iz + 1] +
                       rptr[ixs + iyms + iz + 1];

                // First
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (rptr[ixs + iys + iz - 1] +
                            rptr[ixs + iys + iz + 1] +
                            rz) + c000 * rptr[ixs + iys + iz];


                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c110 * (rptr[ixps + iyps + iz] +
                            rptr[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz] +
                            rzms +
                            rzps);

                // Compute another set of forward edges here
                rzpps = rptr[ixps + iys + iz + 2] +
                        rptr[ixms + iys + iz + 2] +
                        rptr[ixs + iyps + iz + 2] +
                        rptr[ixs + iyms + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (rptr[ixs + iys + iz - 2] +
                            rptr[ixs + iys + iz + 2]) +
                    c200 * (rptr[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + 
                            rptr[ixs + iypps + iz]);


                // Second
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz)] +
                            rptr[ixs + iys + iz + 2] +
                            rzps) + 
                            c000 * rptr[ixs + iys + iz + 1];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    c110 * (rptr[ixps + iyps + iz+1] +
                            rptr[ixps + iyms + iz+1] +
                            rptr[ixms + iyps + iz+1] +
                            rptr[ixms + iyms + iz+1] +
                            rzpps +
                            rz);

                // Compute the forward set of edges here which becomes the trailing set at the top
                rzms = rptr[ixms + iys + iz+3] +
                       rptr[ixps + iys + iz+3] +
                       rptr[ixs + iyms + iz+3] +
                       rptr[ixs + iyps + iz+3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    c200 * (rptr[ixs + iys + iz - 1] +
                            rptr[ixs + iys + iz + 3]) +
                    c200 * (rptr[ixmms + iys + iz+1] +
                            rptr[ixpps + iys + iz+1] +
                            rptr[ixs + iymms + iz+1] + 
                            rptr[ixs + iypps + iz+1]);

                // Third
                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] =
                    c100 * (rptr[ixs + iys + iz+1] +
                            rptr[ixs + iys + iz + 3] +
                            rzpps) +
                            c000 * rptr[ixs + iys + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    c110 * (rptr[ixps + iyps + iz+2] +
                            rptr[ixps + iyms + iz+2] +
                            rptr[ixms + iyps + iz+2] +
                            rptr[ixms + iyms + iz+2] +
                            rzms +
                            rzps);

                // Compute the forward set of edges here which becomes the middle set at the top
                rz = rptr[ixps + iys + iz + 4] +
                     rptr[ixms + iys + iz + 4] +
                     rptr[ixs + iyps + iz + 4] +
                     rptr[ixs + iyms + iz + 4];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    c200 * (rptr[ixs + iys + (iz)] +
                            rptr[ixs + iys + (iz + 4)]) +
                    c200 * (rptr[ixmms + iys + iz+2] +
                            rptr[ixpps + iys + iz+2] +
                            rptr[ixs + iymms + iz+2] + 
                            rptr[ixs + iypps + iz+2]);

                // Fourth
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] =
                    c100 * (rptr[ixs + iys + iz+2] +
                            rptr[ixs + iys + iz + 4] +
                            rzms) + 
                            c000 * rptr[ixs + iys + iz + 3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] +=
                    c110 * (rptr[ixps + iyps + iz+3] +
                            rptr[ixps + iyms + iz+3] +
                            rptr[ixms + iyps + iz+3] +
                            rptr[ixms + iyms + iz+3] +
                            rzpps +
                            rz);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] +=
                    c200 * (rptr[ixs + iys + iz+1] +
                            rptr[ixs + iys + iz + 5]) +
                    c200 * (rptr[ixmms + iys + iz+3] +
                            rptr[ixpps + iys + iz+3] +
                            rptr[ixs + iymms + iz+3] + 
                            rptr[ixs + iypps + iz+3]);


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}


template <typename RmgType>
rmg_double_t FD_app_del2c (RmgType * a, RmgType * b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int ix, iy, iz, ibrav;
    int incy, incx;
    rmg_double_t cc = 0.0, fcx, fcy, fcz, fc, fc1, fc2;
    int ixs, iys, ixms, ixps, iyms, iyps;

    BaseGrid G;

    ibrav = G.get_ibrav_type();

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    switch (ibrav)
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:

        if (G.get_anisotropy() < 1.0000001)
        {

            cc = -2.0 / (gridhx * gridhx * get_xside() * get_xside());
            cc = cc - 2.0 / (gridhy * gridhy * get_yside() * get_yside());
            cc = cc - 2.0 / (gridhz * gridhz * get_zside() * get_zside());
            fcx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());

            for (ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                if (dimy % 2)
                {

                    for (iy = 1; iy <= dimy; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (iz = 1; iz <= dimz; iz++)
                        {

                            b[ixs + iys + iz] =
                                cc * a[ixs + iys + iz] +
                                fcx * (a[ixms + iys + iz] +
                                       a[ixps + iys + iz] +
                                       a[ixs + iyms + iz] +
                                       a[ixs + iyps + iz] +
                                       a[ixs + iys + (iz - 1)] + a[ixs + iys + (iz + 1)]);


                        }       /* end for */

                    }           /* end for */

                }
                else
                {

                    for (iy = 1; iy <= dimy; iy += 2)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (iz = 1; iz <= dimz; iz++)
                        {

                            b[ixs + iys + iz] =
                                cc * a[ixs + iys + iz] +
                                fcx * (a[ixs + iys + (iz - 1)] +
                                       a[ixs + iys + (iz + 1)] +
                                       a[ixms + iys + iz] +
                                       a[ixps + iys + iz] +
                                       a[ixs + iyms + iz] + a[ixs + iyps + iz]);

                            b[ixs + iyps + iz] =
                                cc * a[ixs + iyps + iz] +
                                fcx * (a[ixs + iyps + (iz - 1)] +
                                       a[ixs + iyps + (iz + 1)] +
                                       a[ixms + iyps + iz] +
                                       a[ixps + iyps + iz] +
                                       a[ixs + iys + iz] + a[ixs + iyps + incy + iz]);

                        }       /* end for */

                    }           /* end for */

                }               /* end if */

            }                   /* end for */

        }
        else
        {

            cc = -2.0 / (gridhx * gridhx * get_xside() * get_xside());
            cc = cc - 2.0 / (gridhy * gridhy * get_yside() * get_yside());
            cc = cc - 2.0 / (gridhz * gridhz * get_zside() * get_zside());
            fcx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
            fcy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
            fcz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

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

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fcx * a[ixms + iys + iz] +
                            fcx * a[ixps + iys + iz] +
                            fcy * a[ixs + iyms + iz] +
                            fcy * a[ixs + iyps + iz] +
                            fcz * a[ixs + iys + (iz - 1)] + fcz * a[ixs + iys + (iz + 1)];


                    }           /* end for */

                }               /* end for */

            }                   /* end for */

        }                       /* end if */

        break;

    case CUBIC_BC:

        cc = -2.0 / (gridhx * gridhx * get_xside() * get_xside());
        fc = 1.0 / (4.0 * gridhx * gridhx * get_xside() * get_xside());

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

                    b[ixs + iys + iz] =
                        cc * a[ixs + iys + iz] +
                        fc * a[ixms + iyms + iz - 1] +
                        fc * a[ixms + iys + iz] +
                        fc * a[ixs + iyms + iz] +
                        fc * a[ixs + iys + iz - 1] +
                        fc * a[ixs + iys + iz + 1] +
                        fc * a[ixs + iyps + iz] +
                        fc * a[ixps + iys + iz] + fc * a[ixps + iyps + iz + 1];


                }               /* end for */

            }                   /* end for */

        }                       /* end for */
        break;

    case CUBIC_FC:

        cc = -6.0 / (gridhx * gridhx * get_xside() * get_xside());
        fc = 1.0 / (2.0 * gridhx * gridhx * get_xside() * get_xside());

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

                    b[ixs + iys + iz] =
                        cc * a[ixs + iys + iz] +
                        fc * a[ixms + iys + iz] +
                        fc * a[ixms + iys + iz + 1] +
                        fc * a[ixms + iyps + iz] +
                        fc * a[ixs + iyms + iz] +
                        fc * a[ixs + iyms + iz + 1] +
                        fc * a[ixs + iys + iz - 1] +
                        fc * a[ixs + iys + iz + 1] +
                        fc * a[ixs + iyps + iz - 1] +
                        fc * a[ixs + iyps + iz] +
                        fc * a[ixps + iyms + iz] +
                        fc * a[ixps + iys + iz - 1] + fc * a[ixps + iys + iz];


                }               /* end for */

            }                   /* end for */

        }                       /* end for */
        break;

    case HEXAGONAL:

        cc = -4.0 / (gridhx * gridhx * get_xside() * get_xside());
        cc = cc - 2.0 / (gridhz * gridhz * get_zside() * get_zside());
        fc1 = 2.0 / (3.0 * gridhx * gridhx * get_xside() * get_xside());
        fc2 = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

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

                    b[ixs + iys + iz] =
                        cc * a[ixs + iys + iz] +
                        fc1 * a[ixps + iys + iz] +
                        fc1 * a[ixps + iyms + iz] +
                        fc1 * a[ixs + iyms + iz] +
                        fc1 * a[ixms + iys + iz] +
                        fc1 * a[ixms + iyps + iz] +
                        fc1 * a[ixs + iyps + iz] +
                        fc2 * a[ixs + iys + iz + 1] + fc2 * a[ixs + iys + iz - 1];


                }               /* end for */

            }                   /* end for */

        }                       /* end for */
        break;

    default:
        rmg_error_handler ("Lattice type not implemented");

    }                           /* end switch */


    /* Return the diagonal component of the operator */
    return cc;

}                               /* end app_del2c */


template <typename RmgType>
rmg_double_t FD_app_cil_fourth_standard (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int ibrav;
    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    rmg_double_t ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    rmg_double_t ihx, ihy, ihz;
    BaseGrid G;

    ibrav = G.get_ibrav_type();

    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        rmg_error_handler("Grid symmetry not programmed yet in FD_app_cil_fourth_standard.\n");
    }

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;


    if (G.get_anisotropy() < 1.000001)
    {

        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
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

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        cc * rptr[ixs + iys + iz] +
                        fcx * (rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] +
                                rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        ecxy * (rptr[ixms + iys + iz - 1] +
                                rptr[ixps + iys + iz - 1] +
                                rptr[ixs + iyms + iz - 1] +
                                rptr[ixs + iyps + iz - 1] +
                                rptr[ixms + iyms + iz] +
                                rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] +
                                rptr[ixps + iyps + iz] +
                                rptr[ixms + iys + iz + 1] +
                                rptr[ixps + iys + iz + 1] +
                                rptr[ixs + iyms + iz + 1] + rptr[ixs + iyps + iz + 1]);


                }           /* end for */

            }               /* end for */

        }                   /* end for */
    }
    else
    {

        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
        ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

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

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        cc * rptr[ixs + iys + iz] +
                        fcx * rptr[ixms + iys + iz] +
                        fcx * rptr[ixps + iys + iz] +
                        fcy * rptr[ixs + iyms + iz] +
                        fcy * rptr[ixs + iyps + iz] +
                        fcz * rptr[ixs + iys + (iz - 1)] + fcz * rptr[ixs + iys + (iz + 1)];

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
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
                        ecyz * rptr[ixs + iyms + iz + 1] + ecyz * rptr[ixs + iyps + iz + 1];


                }           /* end for */

            }               /* end for */

        }                   /* end for */

    }                       /* end if */



    return cc;

}                               /* end app_cil */



template <typename RmgType>
rmg_double_t FD_app_cil_fourth_global (RmgType * rptr, RmgType * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t c000, c100, c110;
    rmg_double_t ihx;

    incx = (FIXED_ZDIM + 2) * (FIXED_YDIM + 2);
    incy = FIXED_ZDIM + 2;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    c000 = (-4.0 / 3.0) * (ihx + ihx + ihx);
    c100 = (5.0 / 6.0) * ihx + (c000 / 8.0);
    c110 = (1.0 / 12.0) * (ihx + ihx);


    // Handle the general case first
    for (ix = 1; ix < FIXED_XDIM + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < FIXED_YDIM + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < FIXED_ZDIM + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
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

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    return c000;

}


template <typename RmgType>
void FD_app_cir_fourth_standard (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ibrav;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t c000, c100;
    BaseGrid G;

    ibrav = G.get_ibrav_type();

    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        rmg_error_handler("Grid symmetry not programmed yet in app_cir_fourth.\n");
    }

    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    c000 = 0.5;
    c100 = 1.0 / 12.0;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < dimz + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + 
                    c000 *  rptr[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}

template <typename RmgType>
void FD_app_cir_fourth_global (RmgType * rptr, RmgType * b)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t c000, c100;

    incx = (FIXED_ZDIM + 2) * (FIXED_YDIM + 2);
    incy = FIXED_ZDIM + 2;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    c000 = 0.5;
    c100 = 1.0 / 12.0;


    for (ix = 1; ix < FIXED_XDIM + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < FIXED_YDIM + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < FIXED_ZDIM + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + 
                    c000 * rptr[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}

// Wrappers to call these from C
extern "C" double FD_app_cil_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_sixth_standard<double> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_sixth_standard<float> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}

extern "C" double FD_app_cil_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_sixth_global<double> (rptr, b, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_sixth_global<float> (rptr, b, gridhx, gridhy, gridhz);

}

extern "C" void FD_app_cir_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz)
{

    FD_app_cir_sixth_standard<double> (rptr, b, dimx, dimy, dimz);

}
extern "C" void FD_app_cir_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz)
{

    FD_app_cir_sixth_standard<float> (rptr, b, dimx, dimy, dimz);

}
extern "C" void FD_app_cir_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b)
{

    FD_app_cir_sixth_global<double> (rptr, b);

}
extern "C" void FD_app_cir_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b)
{

    FD_app_cir_sixth_global<float> (rptr, b);

}
extern "C" double app_del2c(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_del2c<double> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double app_del2c_f(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_del2c<float> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_del2c_rmg_complex(complex<double> *rptr, complex<double> *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_del2c<complex <double> > (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_fourth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_fourth_standard<double> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_fourth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_fourth_standard<float> (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_fourth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_fourth_global<double> (rptr, b, gridhx, gridhy, gridhz);

}
extern "C" double FD_app_cil_fourth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_cil_fourth_global<float> (rptr, b, gridhx, gridhy, gridhz);

}
extern "C" void FD_app_cir_fourth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz)
{

    FD_app_cir_fourth_standard<double> (rptr, b, dimx, dimy, dimz);

}
extern "C" void FD_app_cir_fourth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz)
{

    FD_app_cir_fourth_standard<float> (rptr, b, dimx, dimy, dimz);

}
extern "C" void FD_app_cir_fourth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b)
{

    FD_app_cir_fourth_global<double> (rptr, b);

}
extern "C" void FD_app_cir_fourth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b)
{

    FD_app_cir_fourth_global<float> (rptr, b);

}
