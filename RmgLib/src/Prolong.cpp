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

// This is a high order interpolation/prolongation operators used when constructing the
// charge density on the high density grid.
template void Prolong::prolong (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);

Prolong::Prolong(int ratio_in, int order_in, TradeImages &TR_in) : ratio(ratio_in), order(order_in), TR(TR_in)
{
    /*Order has to be even number */
    if (order % 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for even orders.");


    double c[MAX_PROLONG_ORDER];

    for (int ix = 0; ix < MAX_PROLONG_ORDER; ix++)
    {
        for (int iy = 0; iy < MAX_PROLONG_ORDER; iy++)
        {
            a[ix][iy] = 0.0;
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        double fraction = (double) i / (double)ratio;
        cgen_prolong (c, fraction);

        for (int iy = 0; iy < order; iy++)
        {
            int k = iy + (MAX_PROLONG_ORDER - order) / 2;
            a[i][k] = c[iy];
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
    T *sg_half = new T[(half_dimx + 10) * (half_dimy + 10) * (half_dimz + 10)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, 5, FULL_TRADE);

    int incy = dimz / ratio + 10;
    int incx = (dimz / ratio + 10) * (dimy / ratio + 10);

    int incy2 = dimz;
    int incx2 = dimz * dimy;

    int incx3 = (dimz / ratio + 10) * dimy;



    T* fulla = new T[dimx * (dimy / ratio + 10) * (dimz / ratio + 10)];
    T* fullb = new T[ dimx * dimy * (dimz / ratio + 10)];


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx / ratio; ix++)
        {
            for (int iy = 0; iy < dimy / ratio + 10; iy++)
            {
                for (int iz = 0; iz < dimz / ratio + 10; iz++)
                {
                    fulla[((ratio * ix) + i) * incx + iy * incy + iz] =
                        a[i][0] * sg_half[(ix + 1) * incx + iy * incy + iz] +
                        a[i][1] * sg_half[(ix + 2) * incx + iy * incy + iz] +
                        a[i][2] * sg_half[(ix + 3) * incx + iy * incy + iz] +
                        a[i][3] * sg_half[(ix + 4) * incx + iy * incy + iz] +
                        a[i][4] * sg_half[(ix + 5) * incx + iy * incy + iz] +
                        a[i][5] * sg_half[(ix + 6) * incx + iy * incy + iz] +
                        a[i][6] * sg_half[(ix + 7) * incx + iy * incy + iz] +
                        a[i][7] * sg_half[(ix + 8) * incx + iy * incy + iz] +
                        a[i][8] * sg_half[(ix + 9) * incx + iy * incy + iz] +
                        a[i][9] * sg_half[(ix + 10) * incx + iy * incy + iz];
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
                for (int iz = 0; iz < dimz / ratio + 10; iz++)
                {
                    fullb[ix * incx3 + (ratio * iy + i) * incy + iz] =
                        a[i][0] * fulla[ix * incx + (iy + 1) * incy + iz] +
                        a[i][1] * fulla[ix * incx + (iy + 2) * incy + iz] +
                        a[i][2] * fulla[ix * incx + (iy + 3) * incy + iz] +
                        a[i][3] * fulla[ix * incx + (iy + 4) * incy + iz] +
                        a[i][4] * fulla[ix * incx + (iy + 5) * incy + iz] +
                        a[i][5] * fulla[ix * incx + (iy + 6) * incy + iz] +
                        a[i][6] * fulla[ix * incx + (iy + 7) * incy + iz] +
                        a[i][7] * fulla[ix * incx + (iy + 8) * incy + iz] +
                        a[i][8] * fulla[ix * incx + (iy + 9) * incy + iz] +
                        a[i][9] * fulla[ix * incx + (iy + 10) * incy + iz];
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
                    full[ix * incx2 + iy * incy2 + ratio * iz + i] =
                        a[i][0] * fullb[ix * incx3 + iy * incy + iz + 1] +
                        a[i][1] * fullb[ix * incx3 + iy * incy + iz + 2] +
                        a[i][2] * fullb[ix * incx3 + iy * incy + iz + 3] +
                        a[i][3] * fullb[ix * incx3 + iy * incy + iz + 4] +
                        a[i][4] * fullb[ix * incx3 + iy * incy + iz + 5] +
                        a[i][5] * fullb[ix * incx3 + iy * incy + iz + 6] +
                        a[i][6] * fullb[ix * incx3 + iy * incy + iz + 7] +
                        a[i][7] * fullb[ix * incx3 + iy * incy + iz + 8] +
                        a[i][8] * fullb[ix * incx3 + iy * incy + iz + 9] +
                        a[i][9] * fullb[ix * incx3 + iy * incy + iz + 10];
                }
            }
        }
    }


    delete [] fulla;
    delete [] fullb;
    delete [] sg_half;

}

