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
#include "make_conf.h"
#include "common_prototypes.h"
#include "transition.h"

#define TWO 2.0


// Gamma point float version
void CPP_genvpsi (float * psi, float * sg_twovpsi, double * vtot, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{
    for(int ix = 0;ix < dimx*dimy*dimz;ix++) sg_twovpsi[ix] = TWO * psi[ix] * vtot[ix];
}                               /* end genvpsi */

// Gamma point double version
void CPP_genvpsi (double * psi, double * sg_twovpsi, double * vtot, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{
    for(int ix = 0;ix < dimx*dimy*dimz;ix++) sg_twovpsi[ix] = TWO * psi[ix] * vtot[ix];
}

// complex float version
void CPP_genvpsi (std::complex<float> * psi, std::complex<float> * sg_twovpsi, double * vtot, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    std::complex<double> *kd = (std::complex<double> *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[ix * incx + iy * incy + iz] = (std::complex<float>)
                    (
                    2.0 * (std::complex<double>)psi[ix * incx + iy * incy + iz] *
                    (vtot[ix * incx + iy * incy + iz] + 0.5 * kmag) +
                    2.0 * kd[ix * incx + iy * incy + iz]);
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end genvpsi */


// complex double version
void CPP_genvpsi (std::complex<double> * psi, std::complex<double> * sg_twovpsi, double * vtot, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    std::complex<double> *kd = (std::complex<double> *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[ix * incx + iy * incy + iz] = (std::complex<double>)
                    (
                    2.0 * (std::complex<double>)psi[ix * incx + iy * incy + iz] *
                    (vtot[ix * incx + iy * incy + iz] + 0.5 * kmag) +
                    2.0 * kd[ix * incx + iy * incy + iz]);
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end genvpsi */

