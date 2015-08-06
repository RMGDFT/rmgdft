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


#include <stdio.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"

void write_avgd (double * rho)
{

    int ix, iy, iz, poff;
    int px, py, pz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FNZ_GRID;

    double t1;
    double *zvec;
    double hzzgrid;

    hzzgrid = get_hzzgrid();

    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FNZ_GRID = get_FNZ_GRID();

    my_malloc(zvec, FNZ_GRID, double);

    /* Get this processors offset */
    pe2xyz (pct.gridpe, &px, &py, &pz);
    poff = pz * FPZ0_GRID;


    /* Zero out result vector */
    for (iz = 0; iz < FNZ_GRID; iz++)
        zvec[iz] = ZERO;


    /* Loop over this processor */
    for (iz = 0; iz < FPZ0_GRID; iz++)
    {

        t1 = ZERO;
        for (ix = 0; ix < FPX0_GRID; ix++)
        {

            for (iy = 0; iy < FPY0_GRID; iy++)
            {

                t1 += rho[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];

            }                   /* end for */

        }                       /* end for */

        t1 = t1 * get_vel() / hzzgrid;

        zvec[iz + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    iz = FNZ_GRID;
    global_sums (zvec, &iz, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf ("\n\n Planar average of the electrostatic density\n");
        for (iz = 0; iz < FNZ_GRID; iz++)
        {
            t1 = iz * hzzgrid;
            printf (" %f %f\n", t1, zvec[iz]);
        }
        fflush (NULL);
    }

    free(zvec);
}                               /* end get_avgd */

/******/
