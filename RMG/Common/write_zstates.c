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

void write_zstates (STATE * states)
{

rmg_error_handler("Requires updating.");

#if 0
    int ix, iy, iz, poff, PX0_GRID, PY0_GRID, PZ0_GRID;
    int px, py, pz;
    int istate;
    int idx, incix, inciy;
    char newname[MAX_PATH + 20];
    STATE *sp;
    double t1;
    double *zvec;
    double *tmp_psi;
    FILE *avg;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();

    my_malloc (tmp_psi, get_P0_BASIS(), double);
    my_malloc (zvec, get_NZ_GRID(), double);
    /* Get this processors offset */
    pe2xyz (pct.gridpe, &px, &py, &pz);
    poff = pz * PZ0_GRID;

    incix = PY0_GRID * PZ0_GRID;
    inciy = PZ0_GRID;
    for (istate = 0; istate < ct.num_states; istate++)
    {


        sp = &states[istate];

        /* Zero out result vector */
        for (iz = 0; iz < get_NZ_GRID(); iz++)
            zvec[iz] = 0.0;


        gather_psi (tmp_psi, NULL, sp, 0);
        /* Loop over this processor */
        for (iz = 0; iz < PZ0_GRID; iz++)
        {
            t1 = 0.0;
            for (ix = 0; ix < PX0_GRID; ix++)
            {
                for (iy = 0; iy < PY0_GRID; iy++)
                {
                    idx = ix * incix + iy * inciy + iz;
                    t1 += tmp_psi[idx] * tmp_psi[idx];
                }
            }

            zvec[iz + poff] = t1 * get_vel() / (get_hzgrid() * 4.0 * PI);
        }

        /* Now sum over all processors */

        iz = get_NZ_GRID();
        global_sums (zvec, &iz, pct.grid_comm);

        if (pct.gridpe == 0)
        {
            /* Make the new output file name */
            sprintf (newname, "%s%d", "zavg.", istate);
            my_fopen (avg, newname, "w+");

            for (iz = 0; iz < get_NZ_GRID(); iz++)
            {
                t1 = iz * get_hzgrid() * get_zside();
                fprintf (avg, "%f %f\n", t1, zvec[iz]);
            }
            fprintf
                (avg, "\n\n aaa%daaa Planar average of state %d with energy %7.2f eV\n",
                 istate, istate, sp->eig[0] * Ha_eV);
            fflush (NULL);
        }
        fclose (avg);
    }                           /* istate */

    my_free(zvec);
#endif
}                               /* end write_zstates */

/******/
