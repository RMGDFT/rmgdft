/*
 *
 * Copyright (c) 2014, Emil Briggs
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

template double ApplyAOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double ApplyAOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {

        return CPP_app_cil_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);

    }
    else if(ct.discretization == CENTRAL_DISCRETIZATION) {

        double cc = CPP_app_del2_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
        return cc;

    }
    
    rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method."); 
    return false;
}

