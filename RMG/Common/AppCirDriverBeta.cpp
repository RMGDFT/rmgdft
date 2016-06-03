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
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "transition.h"

template void AppCirDriverBeta<double>(Lattice *, TradeImages *, double *, double *, int, int, int, int);
template void AppCirDriverBeta<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, int);


template <typename RmgType>
void AppCirDriverBeta (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order)
{

    RmgTimer RT("App_cir");
    FiniteDiff FD(L);;
    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    RmgType *rptr = new RmgType[sbasis + 64]();

    int incx = dimy * dimz;
    int incy = dimz;

    if(order == APP_CI_FOURTH) {
        int incx1 = (dimy + 2) * (dimz + 2);
        int incy1 = (dimz + 2);
        for(int ix=0;ix < dimx;ix++) {
            for(int iy=0;iy < dimy;iy++) {
                for(int iz=0;iz < dimz;iz++) {
                    rptr[(ix + 1)*incx1 + (iy + 1)*incy1 + iz + 1] = a[ix*incx + iy*incy + iz];
                }
            }
        }

        if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {
            FD.app_cir_fourth (rptr, b, dimx, dimy, dimz);
        }
        else if(ct.discretization == CENTRAL_DISCRETIZATION) {
            for(int ix=0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];
        }
        else {
             rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method.");
        }

    }
    else if(order == APP_CI_SIXTH) {
        int incx1 = (dimy + 4) * (dimz + 4);
        int incy1 = (dimz + 4);
        for(int ix=0;ix < dimx;ix++) {
            for(int iy=0;iy < dimy;iy++) {
                for(int iz=0;iz < dimz;iz++) {
                    rptr[(ix + 2)*incx1 + (iy + 2)*incy1 + iz + 2] = a[ix*incx + iy*incy + iz];
                }
            }
        }

        if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {
            FD.app_cir_sixth (rptr, b, dimx, dimy, dimz);
        }
        else if(ct.discretization == CENTRAL_DISCRETIZATION) {
            for(int ix=0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];
        }
        else {
             rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method.");
        }

    }
    else if(order == APP_CI_EIGHT) {

        if(ct.discretization != CENTRAL_DISCRETIZATION)
            rmg_error_handler(__FILE__, __LINE__, "Eighth order kohn-sham finite differencing is only supported for central fd operators.");

        for(int ix=0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];

    }
    else if(order == APP_CI_TEN) {

        if(ct.discretization != CENTRAL_DISCRETIZATION)
            rmg_error_handler(__FILE__, __LINE__, "Tenth order kohn-sham finite differencing is only supported for central fd operators.");

        for(int ix=0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];

    }
    else if(order == APP_CI_TWELVE) {

        if(ct.discretization != CENTRAL_DISCRETIZATION)
            rmg_error_handler(__FILE__, __LINE__, "Twelfth order kohn-sham finite differencing is only supported for central fd operators.");

        for(int ix=0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];

    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIR order not programmed yet in AppCirDriverBeta.\n");
    }

    delete [] rptr;
    return;

}


