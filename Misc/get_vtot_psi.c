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

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "common_prototypes.h"

void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio)
{
    int idx, ione = 1;
    int ix, iy,iz, idx1;

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        idx = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
        QMD_dcopy (idx, vtot, ione, vtot_psi, ione);
    }
    else if(grid_ratio == 2) {

        FftFilterFine(vtot,  1.0 / (double)grid_ratio );

        for(ix = 0; ix < get_FPX0_GRID()/2; ix++)
        for(iy = 0; iy < get_FPY0_GRID()/2; iy++)
        for(iz = 0; iz < get_FPZ0_GRID()/2; iz++)
        {
            idx = ix * get_FPY0_GRID()/2 * get_FPZ0_GRID()/2 + iy * get_FPZ0_GRID()/2 + iz;
            idx1 = 2 *ix * get_FPY0_GRID() * get_FPZ0_GRID() + 2 *iy * get_FPZ0_GRID() + 2*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
    else if(grid_ratio == 3) {

        FftFilterFine(vtot,  1.0 / (double)grid_ratio );

        for(ix = 0; ix < get_FPX0_GRID()/3; ix++)
        for(iy = 0; iy < get_FPY0_GRID()/3; iy++)
        for(iz = 0; iz < get_FPZ0_GRID()/3; iz++)
        {
            idx = ix * get_FPY0_GRID()/3 * get_FPZ0_GRID()/3 + iy * get_FPZ0_GRID()/3 + iz;
            idx1 = 3 *ix * get_FPY0_GRID() * get_FPZ0_GRID() + 3 *iy * get_FPZ0_GRID() + 3*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
    else
    {
        mg_restrict_6 (vtot, vtot_psi, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), grid_ratio);
    }
}
