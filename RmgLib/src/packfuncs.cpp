/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/



#include "packfuncs.h"
#include <complex>


template void CPP_pack_stop<double>(double*, double*, int, int, int);
template void CPP_pack_stop<float>(float*, float*, int, int, int);
template void CPP_pack_stop<std::complex<float> >(std::complex<float>*, std::complex<float>*, int, int, int);
template void CPP_pack_stop<std::complex<double> >(std::complex<double>*, std::complex<double>*, int, int, int);
template void CPP_pack_ptos<double>(double*, double*, int, int, int);
template void CPP_pack_ptos<float>(float*, float*, int, int, int);
template void CPP_pack_ptos<std::complex<float> >(std::complex<float> *, std::complex<float>*, int, int, int);
template void CPP_pack_ptos<std::complex<double> >(std::complex<double> *, std::complex<double>*, int, int, int);

template <typename RmgType>
void CPP_pack_stop (RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
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
            for(int idx = 0;idx < dimz;idx++) pg[ix * incx + iy * incy + idx] = sg[ixh * incxs + iyh * incys + 1 + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */

void CPP_pack_stop_convert (float * sg, double * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
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
            for(int idx = 0;idx < dimz;idx++) pg[ix * incx + iy * incy + idx] = (double)sg[ixh * incxs + iyh * incys + 1 + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */

template <typename RmgType>
void CPP_pack_ptos(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
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
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */

void CPP_pack_ptos_convert(float * sg, double * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
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
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = (float)pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */

