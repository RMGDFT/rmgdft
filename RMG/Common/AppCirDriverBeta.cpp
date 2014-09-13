/*
 *
 * Copyright (c) 2013, Emil Briggs
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
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIR order not programmed yet in AppCirDriverBeta.\n");
    }

    delete [] rptr;
    return;

}


