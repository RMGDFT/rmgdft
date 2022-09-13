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
#include <alloca.h>
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "LaplacianCoeff.h"
#include "rmg_error.h"

template double CPP_app_del2_driver<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double CPP_app_del2_driver<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template double CPP_app_del2_driver<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int, bool);

template double CPP_app_del2_driver_int<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver_int<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver_int<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int, bool);
template double CPP_app_del2_driver_int<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int, bool);

template <typename RmgType>
double CPP_app_del2_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{
    return CPP_app_del2_driver_int(L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order, false);
}

template <typename RmgType>
double CPP_app_del2_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, bool alt_flag)
{
    return CPP_app_del2_driver_int(L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order, alt_flag);
}

template <typename RmgType>
double CPP_app_del2_driver_int (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order, bool alt_flag)
{

    double cc = 0.0;
    FiniteDiff FD(L, alt_flag);
    int ibrav = L->get_ibrav_type();
    int sbasis = (dimx + order) * (dimy + order) * (dimz + order);
    int images = order / 2;
    size_t alloc = (sbasis + 64) * sizeof(RmgType);
    RmgType *rptr;
    int special = ((ibrav == HEXAGONAL) ||
                   (ibrav == HEXAGONAL2) ||
                   (ibrav == ORTHORHOMBIC_PRIMITIVE) ||
                   (ibrav == CUBIC_PRIMITIVE) ||
                   (ibrav == TETRAGONAL_PRIMITIVE));


    
    // while alloca is dangerous it's very fast for small arrays and the 64k default limit
    // is fine for linux and 64bit power
    if(alloc <= (size_t)FiniteDiff::allocation_limit)
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
        if(a == NULL) return LC->coeff_and_index[0].coeff;
        cc = FiniteDiffLap (rptr, b, dimx, dimy, dimz, LC);
    }
    else if(!special && (dimx == hdim[0]) && (dimy == hdim[1]) && (dimz == hdim[2]))
    {
        if(a == NULL) return HLC->coeff_and_index[0].coeff;
        cc = FiniteDiffLap (rptr, b, dimx, dimy, dimz, HLC);
    }
    else
    {
        if(order == APP_CI_SECOND) {
            cc = FD.app2_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }
        else if(order == APP_CI_EIGHT) {
            cc = FD.app8_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }
        else if(order == APP_CI_TEN) {
            cc = FD.app10_del2 (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }
        else {

            rmg_error_handler (__FILE__, __LINE__, "APP_DEL2 order not programmed yet in app_del2_driver.\n");
            return 0;   // Just to keep the compiler from complaining

        }
    }

    if(alloc > (size_t)FiniteDiff::allocation_limit) delete [] rptr;

    return cc;

}

