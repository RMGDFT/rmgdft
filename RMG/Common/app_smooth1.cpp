
#include "rmgtypes.h"


/**

 Smooths f and returns result in work
*/
template <typename RmgType>
void CPP_app_smooth1 (RmgType * f, RmgType * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;
    int ixs, iys, ixms, ixps, iyms, iyps;

    rmg_double_t scale, ec, fc, crn, cc, temp;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    cc = 6.0;
    fc = 1.0;
    ec = 0.0;
    crn = 0.0;
    scale = 1.0 / 12.0;


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

                temp = cc * f[ixs + iys + iz] +
                    fc * (f[ixms + iys + iz] +
                          f[ixps + iys + iz] +
                          f[ixs + iyms + iz] +
                          f[ixs + iyps + iz] +
                          f[ixs + iys + (iz - 1)] +
                          f[ixs + iys + (iz + 1)]);

                work[ixs + iys + iz] = scale * temp;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end app_smooth1 */



extern "C" void app_smooth1 (rmg_double_t * f, rmg_double_t * work, int dimx, int dimy, int dimz)
{
    CPP_app_smooth1<rmg_double_t> (f, work, dimx, dimy, dimz);
}

extern "C" void app_smooth1_f (rmg_float_t * f, rmg_float_t * work, int dimx, int dimy, int dimz)
{
    CPP_app_smooth1<rmg_float_t> (f, work, dimx, dimy, dimz);
}
