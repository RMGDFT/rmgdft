/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_zstates.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void write_zstates(STATE *states)
 *   Generates the average of each state along the z-axis.
 *   Integrates over the x and y axes.
 * INPUTS
 *   states: points to orbital structure
 * OUTPUT
 *   print out z-average states
 * PARENTS
 *   main.c 
 * CHILDREN
 *   global_sums.c pe2xyz.c
 * SOURCE
 */


#include <stdio.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"

void write_zstates (STATE * states)
{

#if 1
    int ix, iy, iz, poff, PX0_GRID, PY0_GRID, PZ0_GRID;
    int px, py, pz;
    int istate;
    int idx, incix, inciy;
    char newname[MAX_PATH + 20];
    STATE *sp;
    rmg_double_t t1;
    rmg_double_t *zvec;
    rmg_double_t *tmp_psi;
    FILE *avg;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();

    my_malloc (tmp_psi, get_P0_BASIS(), rmg_double_t);
    my_malloc (zvec, get_NZ_GRID(), rmg_double_t);
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

            zvec[iz + poff] = t1 * ct.vel / (get_hzgrid() * 4.0 * PI);
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
#endif

    my_free(zvec);
}                               /* end write_zstates */

/******/
