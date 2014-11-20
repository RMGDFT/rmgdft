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
#include "const.h"
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "transition.h"

template void ApplyBOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, int);
template void ApplyBOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, int);
template void ApplyBOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, int);
template void ApplyBOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, int);

template <typename RmgType>
void ApplyBOperator (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order)
{

    if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {

        CPP_app_cir_driver (L, T, a, b, dimx, dimy, dimz, order);

    }
    else if(ct.discretization == CENTRAL_DISCRETIZATION) {

        for(int ix = 0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];

    }
    else {

        rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method.");

    }

}
