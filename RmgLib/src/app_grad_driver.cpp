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

template void CPP_app_grad_driver<float>(Lattice *, TradeImages *, float *, float *, float *, float *, int, int, int, double, double, double, int);
template void CPP_app_grad_driver<double>(Lattice *, TradeImages *, double *, double *, double *, double *, int, int, int, double, double, double, int);
template void CPP_app_grad_driver<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template void CPP_app_grad_driver<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template <typename RmgType>
void CPP_app_grad_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * bx, RmgType * by, RmgType * bz, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    RmgTimer RT("App_gradient");
    int sbasis;
    FiniteDiff FD(L);
    sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    RmgType *rptr = new RmgType[sbasis + 64];

    if(order == APP_CI_FOURTH) {

        RmgTimer *RT1 = new RmgTimer("App_gradient: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 2, CENTRAL_TRADE);
        delete(RT1);
        RT1 = new RmgTimer("App_gradient: computation");
        FD.app_gradient_fourth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete(RT1);

    }
    else if(order == APP_CI_SIXTH) {

        RmgTimer *RT1 = new RmgTimer("App_gradient: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 3, CENTRAL_TRADE);
        delete(RT1);
        RT1 = new RmgTimer("App_gradient: computation");
        FD.app_gradient_sixth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete(RT1);

    }
    else if(order == APP_CI_EIGHT) {

        RmgTimer *RT1 = new RmgTimer("App_gradient: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 4, CENTRAL_TRADE);
        delete(RT1);
        RT1 = new RmgTimer("App_gradient: computation");
        FD.app_gradient_eighth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete(RT1);

    }
    else if(order == APP_CI_TEN) {

        RmgTimer *RT1 = new RmgTimer("App_gradient: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 5, CENTRAL_TRADE);
        delete(RT1);
        RT1 = new RmgTimer("App_gradient: computation");
        FD.app_gradient_tenth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete(RT1);

    }
    else if(order == APP_CI_TWELVE) {

        RmgTimer *RT1 = new RmgTimer("App_gradient: trade images");
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, 6, CENTRAL_TRADE);
        delete(RT1);
        RT1 = new RmgTimer("App_gradient: computation");
        FD.app_gradient_twelfth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        delete(RT1);

    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "Finite difference order not programmed yet in app_gradient_driver.\n");
    }



    delete [] rptr;

}

