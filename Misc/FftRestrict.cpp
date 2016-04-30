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

// This function is used to restrict an array from a fine
// grid to a coarser grid.
// If the grids are not identical then frequencies too high to
// to be represented on the coarse grid are fourier filtered before
// a restriction operator is applied.

void FftRestrict (double * fine, double * coarse, int grid_ratio)
{
    int idx;
    int ix, iy,iz, idx1;
 
    int dimx = Rmg_G->get_PX0_GRID(grid_ratio);
    int dimy = Rmg_G->get_PY0_GRID(grid_ratio);
    int dimz = Rmg_G->get_PZ0_GRID(grid_ratio);

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        for(idx = 0;idx < dimx*dimy*dimz;idx++) coarse[idx] = fine[idx];
        return;
    }

    double *temp = new double[dimx*dimy*dimz];
    for(int idx=0;idx<dimx*dimy*dimz;idx++)temp[idx] = fine[idx];

    if(grid_ratio == 2) {

        FftFilter(temp, *fine_pwaves, 1.0 / (double)grid_ratio, LOW_PASS);
        for(ix = 0; ix < dimx/2; ix++)
        for(iy = 0; iy < dimy/2; iy++)
        for(iz = 0; iz < dimz/2; iz++)
        {
            idx = ix * dimy/2 * dimz/2 + iy * dimz/2 + iz;
            idx1 = 2 *ix * dimy * dimz + 2 *iy * dimz + 2*iz;
            coarse[idx] = temp[idx1];
        }

    }
    else if(grid_ratio == 3) {

        FftFilter(temp, *fine_pwaves, 1.0 / (double)grid_ratio, LOW_PASS);
        for(ix = 0; ix < dimx/3; ix++)
        for(iy = 0; iy < dimy/3; iy++)
        for(iz = 0; iz < dimz/3; iz++)
        {
            idx = ix * dimy/3 * dimz/3 + iy * dimz/3 + iz;
            idx1 = 3 *ix * dimy * dimz + 3 *iy * dimz + 3*iz;
            coarse[idx] = temp[idx1];
        }

    }
    else
    {
        FftFilter(temp, *fine_pwaves, 1.0 / (double)grid_ratio, LOW_PASS);
        mg_restrict_6 (temp, coarse, dimx, dimy, dimz, grid_ratio);
    }
    delete [] temp;
}
