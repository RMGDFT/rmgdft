/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "Prolong.h"
#include "blas.h"
#include <complex>
#include <cmath>

// This is a high order interpolation/prolongation operators used when constructing the
// charge density on the high density grid.
template void Prolong::prolong (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);

Prolong::Prolong(int ratio_in, int order_in, TradeImages &TR_in) : ratio(ratio_in), order(order_in), TR(TR_in)
{
    /*Order has to be even number */
    if (order % 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for even orders.");



    for (int ix = 0; ix < MAX_PROLONG_ORDER; ix++)
    {
        for (int iy = 0; iy < MAX_PROLONG_ORDER; iy++)
        {
            a[ix][iy] = 0.0;
            af[ix][iy] = 0.0;
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        double fraction = (double) i / (double)ratio;
        double c[MAX_PROLONG_ORDER] = {0.0};
        cgen_prolong (c, fraction);

        for (int iy = 0; iy < MAX_PROLONG_ORDER; iy++)
        {
            a[i][iy] = c[iy];
            if(fabs(c[iy]) > 1.0e-20) af[i][iy] = (float)c[iy];
        }
    }
}

Prolong::~Prolong(void)
{

}


void Prolong::cgen_prolong (double *coef, double fraction)
{

    double *A = new double[order * order]();
    double *b = new double[order]();
    double *d = new double[order]();
    int *ipvt = new int[order]();
    int info;
    int ione = 1;

    /*filling A , b and d (d is distance from different coarse grids to the interpolated pt, 
       coarse grid spacing is normalized to be ONE) */
    b[0] = 1.0;

    for (int iy = 0; iy < order; iy++)
    {
        d[iy] = iy + 1.0 - fraction - (double) order / 2;
        for (int ix = 0; ix < order; ix++)
        {
            A[iy * order + ix] = pow (d[iy], ix);
        }                       /* end for */
    }                           /* end for */

    /*  solving Ac=b for c using  b = A^(-1) * b  */
    dgesv (&order, &ione, A, &order, ipvt, b, &order, &info);

    for (int ix = 0; ix < order; ix++) coef[ix] = b[ix];

    delete [] A;
    delete [] ipvt;
    delete [] d;
    delete [] b;

}




template <typename T> void Prolong::prolong (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz)
{

    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order;
    int incx = (dimz / ratio + order) * (dimy / ratio + order);

    int incy2 = dimz;
    int incx2 = dimz * dimy;

    int incx3 = (dimz / ratio + order) * dimy;



    T* fulla = new T[dimx * (dimy / ratio + order) * (dimz / ratio + order)];
    T* fullb = new T[ dimx * dimy * (dimz / ratio + order)];

    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx / ratio; ix++)
        {
            for (int iy = 0; iy < dimy / ratio + order; iy++)
            {
                for (int iz = 0; iz < dimz / ratio + order; iz++)
                {
                    T sum = 0.0;
                    T *halfptr = &sg_half[(ix + 1) * incx + iy * incy + iz];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * halfptr[k*incx];
                    fulla[((ratio * ix) + i) * incx + iy * incy + iz] = sum;
                }
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy / ratio; iy++)
            {
                for (int iz = 0; iz < dimz / ratio + order; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fulla[ix * incx + (iy + 1) * incy + iz];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k*incy];
                    fullb[ix * incx3 + (ratio * iy + i) * incy + iz] = sum;
                }
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy; iy++)
            {
                for (int iz = 0; iz < dimz / ratio; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fullb[ix * incx3 + iy * incy + iz + 1];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k];
                    full[ix * incx2 + iy * incy2 + ratio * iz + i] = sum;
                }
            }
        }
    }


    delete [] fulla;
    delete [] fullb;
    delete [] sg_half;

}

