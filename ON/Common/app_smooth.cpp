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



#include <complex>

template <typename RmgType>
void CPP_app_smooth (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

template void CPP_app_smooth<float> (float *, float *, int , int , int );
template void CPP_app_smooth<double> (double *, double *, int , int , int );
template void CPP_app_smooth<std::complex<float> > (std::complex<float> *, std::complex<float> *, int , int , int );
template void CPP_app_smooth<std::complex<double> > (std::complex<double> *, std::complex<double> *, int , int , int );



template <typename RmgType>
void CPP_app_smooth (RmgType * f, RmgType * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;
    int ixs, iys, ixms, ixps, iyms, iyps;

    RmgType scale, ec, fc, crn, cc, temp;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;

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


                temp +=
                    ec * (f[ixms + iys + iz - 1] +
                          f[ixps + iys + iz - 1] +
                          f[ixs + iyms + iz - 1] +
                          f[ixs + iyps + iz - 1] +
                          f[ixms + iyms + iz] +
                          f[ixms + iyps + iz] +
                          f[ixps + iyms + iz] +
                          f[ixps + iyps + iz] +
                          f[ixms + iys + iz + 1] +
                          f[ixps + iys + iz + 1] +
                          f[ixs + iyms + iz + 1] + 
                          f[ixs + iyps + iz + 1]);


                temp +=
                    crn * (f[ixms + iyms + iz - 1] +
                           f[ixms + iyms + iz + 1] +
                           f[ixms + iyps + iz - 1] +
                           f[ixms + iyps + iz + 1] +
                           f[ixps + iyms + iz - 1] +
                           f[ixps + iyms + iz + 1] +
                           f[ixps + iyps + iz - 1] +
                           f[ixps + iyps + iz + 1]);

                work[ixs + iys + iz] = scale * temp;


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}


extern "C" void app_smooth (double * f, double * work, int dimx, int dimy, int dimz)
{
    CPP_app_smooth<double> (f, work, dimx, dimy, dimz);
}

extern "C" void app_smooth_f (float * f, float * work, int dimx, int dimy, int dimz)
{
    CPP_app_smooth<float> (f, work, dimx, dimy, dimz);
}
