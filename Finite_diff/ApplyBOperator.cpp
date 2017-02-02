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
#include "RmgException.h"
#include "rmg_error.h"
#include "transition.h"

// Applies Mehrstellen right hand side operator to a and returns result in b
// The first set of functions takes the input and output grids and a char string that defines
// the grid. More detailed parameters are then passed to the second set which may be accessed
// directly if more control is required.
//
// IN:    Input array a defined on coarse or fine grid
// OUT:   Output array b defined on coarse or fine grid
// IN:    grid = "Coarse" or "Fine" for grid type

template void ApplyBOperator<float>(float *, float *, char *grid);
template void ApplyBOperator<double>(double *, double *, char *grid);
template void ApplyBOperator<std::complex<float> >(std::complex<float> *, std::complex<float> *, char *grid);
template void ApplyBOperator<std::complex<double> >(std::complex<double> *, std::complex<double> *, char *grid);

template void ApplyBOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, int);
template void ApplyBOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, int);
template void ApplyBOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, int);
template void ApplyBOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, int);


template <typename RmgType>
void ApplyBOperator (RmgType * a, RmgType * b, char *grid)
{
    int density;
    const char *coarse = "Coarse";
    const char *fine = "Fine";

    if(!strcmp(grid, coarse)) {
        density = 1;
    }
    else if(!strcmp(grid, fine)) {
        density = Rmg_G->default_FG_RATIO;
    }
    else {
        throw RmgFatalException() << "Error! Grid type " << grid << " not defined in "
                                 << __FILE__ << " at line " << __LINE__ << "\n";
    }

    int dimx = Rmg_G->get_PX0_GRID(density);
    int dimy = Rmg_G->get_PY0_GRID(density);
    int dimz = Rmg_G->get_PZ0_GRID(density);

    ApplyBOperator(&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, ct.kohn_sham_fd_order);

}

template <typename RmgType>
void ApplyBOperator (Lattice *L, TradeImages *T, RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz, int order)
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
