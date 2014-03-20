#include "make_conf.h"
#include "rmgtypes.h"
#include "auxiliary.h"

using namespace std;

template void CPP_pack_stop_axpy<double>(double*, double*, double, int, int, int);
template void CPP_pack_stop_axpy<float>(float*, float*, double, int, int, int);

template <typename RmgType>
void CPP_pack_stop_axpy (RmgType * sg, RmgType * pg, rmg_double_t alpha, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;

    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for (iz = 0; iz < dimz; iz++)
            {

                pg[ix * incx + iy * incy + iz] += alpha * sg[ixh * incxs + iyh * incys + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


} // end CPP_pack_stop_axpy


extern "C" void pack_stop_axpy (rmg_double_t * sg, rmg_double_t * pg, rmg_double_t alpha, int dimx, int dimy, int dimz)
{
    CPP_pack_stop_axpy<rmg_double_t> (sg, pg, alpha, dimx, dimy, dimz);
}


extern "C" void pack_stop_axpy_f (rmg_float_t * sg, rmg_float_t * pg, rmg_double_t alpha, int dimx, int dimy, int dimz)
{
    CPP_pack_stop_axpy<rmg_float_t> (sg, pg, alpha, dimx, dimy, dimz);
}



