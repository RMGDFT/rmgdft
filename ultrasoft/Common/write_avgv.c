/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_avgv.c *****
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
 *   void write_avgv(REAL *vh, REAL *vnuc)
 *   Generates the avergage electrostatic potential along the z-axis.
 *   Integrates over the x and y axes.
 * INPUTS
 *   vh:  Hartree potential
 *   vnuc: ionic potential (pseudopential)
 * OUTPUT
 *   the average is printed out
 * PARENTS
 *   main.c
 * CHILDREN
 *   global_sums.c pe2xyz.c
 * SOURCE
 */



#include <stdio.h>
#include "main.h"


void write_avgv (REAL * vh, REAL * vnuc)
{

    int ix, iy, iz, poff;
    int px, py, pz;
    REAL t1;
    REAL zvec[NZ_GRID];


    /* Get this processors offset */
    pe2xyz (pct.gridpe, &px, &py, &pz);
    poff = pz * PZ0_GRID;


    /* Zero out result vector */
    for (iz = 0; iz < NZ_GRID; iz++)
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

        t1 = t1 * ct.vel / ct.hzgrid / 4.0 / PI;

        zvec[iz + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    iz = NZ_GRID;
    global_sums (zvec, &iz, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf ("\n\n Planar average of the electron potential\n");
        for (iz = 0; iz < NZ_GRID; iz++)
        {
            t1 = iz * ct.hzgrid;
            printf (" %f %f\n", t1, zvec[iz]);
        }

        fflush (NULL);
    }
}                               /* end get_avgv */

/******/
