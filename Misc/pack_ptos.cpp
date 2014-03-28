
#include "auxiliary.h"
#include "BlasWrappers.h"


template <typename RmgType>
void CPP_pack_ptos(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    int ione = 1;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            QMD_copy(dimz, &pg[ix * incx + iy * incy], ione,
                      &sg[ixh * incxs + iyh * incys + 1], ione);

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */


extern "C" void pack_ptos(double * sg, double * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_ptos<double> (sg, pg, dimx, dimy, dimz);
}

extern "C" void pack_ptos_f(float * sg, float * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_ptos<float> (sg, pg, dimx, dimy, dimz);
}
