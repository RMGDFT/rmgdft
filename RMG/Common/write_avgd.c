/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_avgd.c *****
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
 *   void write_avgd(rmg_double_t *rho)
 *   Generates the average electron density along the z-axis.
 *   Integrates over the x and y axes.
 * INPUTS
 *   rho: total valence charge density
 * OUTPUT
 *   the average over xy plane is printed out
 * PARENTS
 *   main.c
 * CHILDREN
 *   globlal_sums.c pe2xyz.c
 * SOURCE
 */


#include <stdio.h>
#include "main.h"

void write_avgd (rmg_double_t * rho)
{

    int ix, iy, iz, poff;
    int px, py, pz;
    rmg_double_t t1;
    rmg_double_t zvec[FNZ_GRID];

    /* Get this processors offset */
    pe2xyz (pct.gridpe, &px, &py, &pz);
    poff = pz * pct.FPZ0_GRID;


    /* Zero out result vector */
    for (iz = 0; iz < FNZ_GRID; iz++)
        zvec[iz] = ZERO;


    /* Loop over this processor */
    for (iz = 0; iz < pct.FPZ0_GRID; iz++)
    {

        t1 = ZERO;
        for (ix = 0; ix < pct.FPX0_GRID; ix++)
        {

            for (iy = 0; iy < pct.FPY0_GRID; iy++)
            {

                t1 += rho[ix * pct.FPY0_GRID * pct.FPZ0_GRID + iy * pct.FPZ0_GRID + iz];

            }                   /* end for */

        }                       /* end for */

        t1 = t1 * ct.vel / ct.hzzgrid;

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
            t1 = iz * ct.hzzgrid;
            printf (" %f %f\n", t1, zvec[iz]);
        }
        fflush (NULL);
    }
}                               /* end get_avgd */

/******/
