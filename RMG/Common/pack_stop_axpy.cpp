/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include "make_conf.h"
#include "packfuncs.h"
#include <complex>
#include <typeinfo>



template void CPP_pack_stop_axpy<double>(double*, double*, double, int, int, int);
template void CPP_pack_stop_axpy<float>(float*, float*, double, int, int, int);
template void CPP_pack_stop_axpy<std::complex<float> >(std::complex<float>*, std::complex<float>*, double, int, int, int);
template void CPP_pack_stop_axpy<std::complex<double> >(std::complex<double>*, std::complex<double>*, double, int, int, int);

template <typename RmgType>
void CPP_pack_stop_axpy (RmgType * sg, RmgType * pg, double alpha, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;

    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);
    RmgType alpha1(alpha);

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for (iz = 0; iz < dimz; iz++)
            {

                pg[ix * incx + iy * incy + iz] = pg[ix * incx + iy * incy + iz] + alpha1 * sg[ixh * incxs + iyh * incys + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


} // end CPP_pack_stop_axpy


extern "C" void pack_stop_axpy (double * sg, double * pg, double alpha, int dimx, int dimy, int dimz)
{
    CPP_pack_stop_axpy<double> (sg, pg, alpha, dimx, dimy, dimz);
}


extern "C" void pack_stop_axpy_f (float * sg, float * pg, double alpha, int dimx, int dimy, int dimz)
{
    CPP_pack_stop_axpy<float> (sg, pg, alpha, dimx, dimy, dimz);
}



