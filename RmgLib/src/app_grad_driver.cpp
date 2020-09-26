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
#include "LaplacianCoeff.h"
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
    RmgType *rptr;
    sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    int images = order / 2;
    int ibrav = L->get_ibrav_type();
    size_t alloc = (sbasis + 64) * sizeof(RmgType);
    int special = ((ibrav == HEXAGONAL) ||
                   (ibrav == HEXAGONAL2) ||
                   (ibrav == ORTHORHOMBIC_PRIMITIVE) || 
                   (ibrav == CUBIC_PRIMITIVE) ||
                   (ibrav == TETRAGONAL_PRIMITIVE));


    // while alloca is dangerous it's very fast for small arrays and the 110k limit
    // is fine for linux and 64bit power
    if(alloc <= 110592)
    {
        rptr = (RmgType *)alloca(alloc);
    }
    else
    {
        rptr = new RmgType[sbasis + 64];
    }

    if(!special || (ibrav == HEXAGONAL) || (ibrav == HEXAGONAL2))
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, FULL_TRADE);
    else
        T->trade_imagesx (a, rptr, dimx, dimy, dimz, images, CENTRAL_TRADE);

    int dim[3]={0,0,0}, hdim[3]={0,0,0};
    if(!special)
    {
        LC->GetDim(dim);
        HLC->GetDim(hdim);
    }
    if(!special && (dimx == dim[0]) && (dimy == dim[1]) && (dimz == dim[2]))
    {
        FiniteDiffGrad (rptr, bx, by, bz, dimx, dimy, dimz, LC);
    }
    else if(!special && (dimx == hdim[0]) && (dimy == hdim[1]) && (dimz == hdim[2]))
    {
        FiniteDiffGrad (rptr, bx, by, bz, dimx, dimy, dimz, HLC);
    }
    else
    {

        if(order == APP_CI_SIXTH) {

            FD.app_gradient_sixth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);

        }
        else if(order == APP_CI_EIGHT) {

            FD.app_gradient_eighth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);

        }
        else if(order == APP_CI_TEN) {

            FD.app_gradient_tenth (rptr, bx, by, bz, dimx, dimy, dimz, gridhx, gridhy, gridhz);

        }
        else {
            rmg_error_handler (__FILE__, __LINE__, "Finite difference order not programmed yet in app_gradient_driver.\n");
        }

    }

    if(alloc > 110592) delete [] rptr;

}

