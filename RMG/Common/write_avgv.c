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


void write_avgv (double * vh, double * vnuc)
{

    int ix, iy, iz, poff;
    int px, py, pz;
    double t1;
    double *zvec;
    int PX0_GRID, PY0_GRID, PZ0_GRID;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();

    my_malloc (zvec, get_NZ_GRID(), double);


    /* Get this processors offset */
    pe2xyz (pct.gridpe, &px, &py, &pz);
    poff = pz * PZ0_GRID;


    /* Zero out result vector */
    for (iz = 0; iz < get_NZ_GRID(); iz++)
        zvec[iz] = ZERO;


    /* Loop over this processor */
    for (iz = 0; iz < PZ0_GRID; iz++)
    {

        t1 = ZERO;
        for (ix = 0; ix < PX0_GRID; ix++)
        {

            for (iy = 0; iy < PY0_GRID; iy++)
            {

                t1 += vh[ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz] -
                    vnuc[ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz];

            }                   /* end for */

        }                       /* end for */

        t1 = t1 * get_vel() / get_hzgrid() / 4.0 / PI;

        zvec[iz + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    iz = get_NZ_GRID();
    global_sums (zvec, &iz, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf ("\n\n Planar average of the electron potential\n");
        for (iz = 0; iz < get_NZ_GRID(); iz++)
        {
            t1 = iz * get_hzgrid();
            printf (" %f %f\n", t1, zvec[iz]);
        }

        fflush (NULL);
    }

    my_free(zvec);
}                               /* end get_avgv */

/******/
