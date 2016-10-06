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
#include <alloca.h>

template double CPP_app_del2_driver<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template <typename RmgType>
double CPP_app_del2_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    RmgTimer RT("App_del2");
    int sbasis;
    RmgType *rptr;
    double cc = 0.0;
    FiniteDiff FD(L);

    if(order == APP_CI_SECOND) {
        sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);
        rptr = (RmgType *)alloca((sbasis + 64) * sizeof(RmgType));
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 1, CENTRAL_TRADE);
        cc = FD.app2_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        //free(rptr);
    }
    else if(order == APP_CI_FOURTH) {
        sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
        rptr = (RmgType *)alloca((sbasis + 64) * sizeof(RmgType));
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 2, CENTRAL_TRADE);
        cc = FD.app4_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        //free(rptr);
    }
    else if(order == APP_CI_SIXTH) {
        sbasis = (dimx + 6) * (dimy + 6) * (dimz + 6);
        rptr = (RmgType *)alloca((sbasis + 64) * sizeof(RmgType));
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 3, CENTRAL_TRADE);
        cc = FD.app6_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        //free(rptr);
    }
    else if(order == APP_CI_EIGHT) {
        sbasis = (dimx + 8) * (dimy + 8) * (dimz + 8);
        rptr = new RmgType [sbasis + 64];
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 4, CENTRAL_TRADE);
        cc = FD.app8_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete [] rptr;
        //free(rptr);
    }
    else if(order == APP_CI_TEN) {
        sbasis = (dimx + 10) * (dimy + 10) * (dimz + 10);
        rptr = (RmgType *)alloca((sbasis + 64) * sizeof(RmgType));
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 5, CENTRAL_TRADE);
        cc = FD.app10_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        //free(rptr);
    }
    else if(order == APP_CI_TWELVE) {
        sbasis = (dimx + 12) * (dimy + 12) * (dimz + 12);
        rptr = (RmgType *)alloca((sbasis + 64) * sizeof(RmgType));
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 6, CENTRAL_TRADE);
        cc = FD.app12_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        //free(rptr);
    }
    else {

        rmg_error_handler (__FILE__, __LINE__, "APP_DEL2 order not programmed yet in app_del2_driver.\n");
        return 0;   // Just to keep the compiler from complaining

    }

    return cc;

}

