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

#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include <complex>

template double CPP_app_cil_driver_threaded<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double CPP_app_cil_driver_threaded<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double CPP_app_cil_driver_threaded<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double CPP_app_cil_driver_threaded<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template <typename RmgType>
double CPP_app_cil_driver_threaded (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    RmgTimer RT("App_cil");
    int sbasis;
    double cc = 0.0;
    FiniteDiff FD(L);
    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    RmgType *rptr = new RmgType[sbasis + 64];

    if(order == APP_CI_FOURTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_fourth_threaded (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_SIXTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_sixth (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIL order not programmed yet in app_cil_driver_threaded.\n");
    }

    delete [] rptr;
    return cc;

}

