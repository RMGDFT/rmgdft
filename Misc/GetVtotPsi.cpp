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
#include <math.h>
#include <float.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "RmgParallelFft.h"

// This function is used to transfer the potential from the fine
// grid to the coarse grid that the wavefunctions are defined on.
// If the grids are not identical then frequencies too high to
// to be represented on the coarse grid are fourier filtered before
// a restriction operator is applied.

void GetVtotPsi (double * vtot_psi, double * in_vtot, int grid_ratio)
{
    int idx;
    int ix, iy,iz, idx1;
 
    int dimx = Rmg_G->get_PX0_GRID(grid_ratio);
    int dimy = Rmg_G->get_PY0_GRID(grid_ratio);
    int dimz = Rmg_G->get_PZ0_GRID(grid_ratio);

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        for(idx = 0;idx < dimx*dimy*dimz;idx++) vtot_psi[idx] = in_vtot[idx];
        return;
    }

    double *vtot = new double[dimx*dimy*dimz];
    for(int idx=0;idx<dimx*dimy*dimz;idx++)vtot[idx] = in_vtot[idx];

    if(grid_ratio == 2) {

        FftFilter(vtot, *fine_pwaves, sqrt(ct.filter_factor) / (double)grid_ratio, LOW_PASS);

        for(ix = 0; ix < dimx/2; ix++)
        for(iy = 0; iy < dimy/2; iy++)
        for(iz = 0; iz < dimz/2; iz++)
        {
            idx = ix * dimy/2 * dimz/2 + iy * dimz/2 + iz;
            idx1 = 2 *ix * dimy * dimz + 2 *iy * dimz + 2*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
    else if(grid_ratio == 3) {

        FftFilter(vtot, *fine_pwaves, sqrt(ct.filter_factor) / (double)grid_ratio, LOW_PASS);

        for(ix = 0; ix < dimx/3; ix++)
        for(iy = 0; iy < dimy/3; iy++)
        for(iz = 0; iz < dimz/3; iz++)
        {
            idx = ix * dimy/3 * dimz/3 + iy * dimz/3 + iz;
            idx1 = 3 *ix * dimy * dimz + 3 *iy * dimz + 3*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
    else if(grid_ratio == 4) {

        FftFilter(vtot, *fine_pwaves, sqrt(ct.filter_factor) / (double)grid_ratio, LOW_PASS);

        for(ix = 0; ix < dimx/4; ix++)
        for(iy = 0; iy < dimy/4; iy++)
        for(iz = 0; iz < dimz/4; iz++)
        {
            idx = ix * dimy/4 * dimz/4 + iy * dimz/4 + iz;
            idx1 = 4 *ix * dimy * dimz + 4 *iy * dimz + 4*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
    else
    {

        FftFilter(vtot, *fine_pwaves, sqrt(ct.filter_factor) / (double)grid_ratio, LOW_PASS);

        mg_restrict_6 (vtot, vtot_psi, dimx, dimy, dimz, grid_ratio);
    }
    delete [] vtot;
}
