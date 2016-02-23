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
#include <string>
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "RmgException.h"
#include "transition.h"

// Applies Mehrstellen right hand side operator to a and returns result in b
// IN:    Input array a defined on coarse or fine grid
// OUT:   Output array b defined on coarse of fine grid
// IN:    grid = "Coarse" or "Fine" for grid type

template <typename RmgType>
void AppCir (RmgType * a, RmgType * b, char * grid);

template void AppCir<float>(float *, float *, char *);
template void AppCir<double>(double *, double *, char *);
template void AppCir<std::complex<float> >(std::complex<float> *, std::complex<float> *, char *);
template void AppCir<std::complex<double> >(std::complex<double> *, std::complex<double> *, char *);


template <typename RmgType>
void AppCir (RmgType * a, RmgType * b, char * grid)
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

    CPP_app_cir_driver (&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, ct.kohn_sham_fd_order);


}

