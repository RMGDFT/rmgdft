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
#include <TradeImages.h>
#include <transition.h>

template <typename RmgType>
void CPP_app_smooth_test (RmgType * f, RmgType *b, int dimx, int dimy, int dimz);

//template void CPP_app_smooth_test<float> (float *, float *, int , int , int );
template void CPP_app_smooth_test<double> (double *, double *, int , int , int );
//template void CPP_app_smooth_test<std::complex<float> > (std::complex<float> *, std::complex<float> *, int , int , int );
template void CPP_app_smooth_test<std::complex<double> > (std::complex<double> *, std::complex<double> *, int , int , int );



template <typename RmgType>
void CPP_app_smooth_test (RmgType * a, RmgType *b, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;

    RmgType scale, ec, fc, crn, cc, temp;
    int sbasis = (dimx + 8) * (dimy + 8) * (dimz + 8);
    RmgType *f = new RmgType[sbasis + 64]();

    incy = dimz + 8;
    incx = (dimz + 8) * (dimy + 8);

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;

    Rmg_T->trade_imagesx (a, f, dimx, dimy, dimz, 4, FULL_TRADE);

    for (ix = 4; ix < dimx + 4; ix++)
    {

        int ixs = ix * incx;
        int ixms = (ix - 1) * incx;
        int ixps = (ix + 1) * incx;
        int ixmms = (ix - 2) * incx;
        int ixpps = (ix + 2) * incx;

        for (iy = 4; iy < dimy + 4; iy++)
        {

            int iys = iy * incy;
            int iyms = (iy - 1) * incy;
            int iyps = (iy + 1) * incy;
            int iymms = (iy - 2) * incy;
            int iypps = (iy + 2) * incy;

            for (iz = 4; iz < dimz + 4; iz++)
            {
#if 1
b[(ix - 4)*dimy*dimz + (iy - 4)*dimz + iz - 4] = 0.0;
for(int ixx = ix-4;ixx <= ix+4;ixx+=4){
  for(int iyy = iy-4;iyy <= iy+4;iyy+=4){
    for(int izz = iz-4;izz <= iz+4;izz+=4){
        double r = (double)((ixx-ix) * (ixx-ix) + (iyy-iy) * (iyy-iy) + (izz-iz) * (izz-iz));
        if(r > 0.1) {
            r = sqrt(r);
            b[(ix - 4)*dimy*dimz + (iy - 4)*dimz + iz - 4] += f[ixx*incx + iyy*incy + izz]/r;
        }
        else {
            b[(ix - 4)*dimy*dimz + (iy - 4)*dimz + iz - 4] += 16.0*f[ixx*incx + iyy*incy + izz];
        }
    }
  }
}
#else

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


                temp +=
                    fc * (f[ixmms + iys + iz] +
                          f[ixpps + iys + iz] +
                          f[ixs + iymms + iz] +
                          f[ixs + iypps + iz] +
                          f[ixs + iys + (iz - 2)] + 
                          f[ixs + iys + (iz + 2)]);

                temp +=
                    ec * (f[ixmms + iys + iz - 2] +
                          f[ixpps + iys + iz - 2] +
                          f[ixs + iymms + iz - 2] +
                          f[ixs + iypps + iz - 2] +
                          f[ixmms + iymms + iz] +
                          f[ixmms + iypps + iz] +
                          f[ixpps + iymms + iz] +
                          f[ixpps + iypps + iz] +
                          f[ixmms + iys + iz + 2] +
                          f[ixpps + iys + iz + 2] +
                          f[ixs + iymms + iz + 2] + 
                          f[ixs + iypps + iz + 2]);

                temp +=
                    crn * (f[ixmms + iymms + iz - 2] +
                           f[ixmms + iymms + iz + 2] +
                           f[ixmms + iypps + iz - 2] +
                           f[ixmms + iypps + iz + 2] +
                           f[ixpps + iymms + iz - 2] +
                           f[ixpps + iymms + iz + 2] +
                           f[ixpps + iypps + iz - 2] +
                           f[ixpps + iypps + iz + 2]);

                b[(ix - 2)*dimy*dimz + (iy - 2)*dimz + iz - 2] = scale * temp;

#endif

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    delete [] f;
}

