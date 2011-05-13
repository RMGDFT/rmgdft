/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <stdio.h>
#include "md.h"

void write_rho_x(REAL * rho, char *ab)
{

    int ix, iy, iz, poff;
    int px, py, pz;
    REAL t1;
    REAL *zvec;
    int pyoff, pzoff;

    my_malloc_init( zvec, FNX_GRID, REAL );
    /* Get this processors offset */
    pe2xyz(pct.gridpe, &px, &py, &pz);
    poff = px * FPX0_GRID;
    pyoff = py * FPY0_GRID;
    pzoff = pz * FPZ0_GRID;



    /* Zero out result vector */
    for (ix = 0; ix < FNX_GRID; ix++)
        zvec[ix] = ZERO;


    /* Loop over this processor */
    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        t1 = 0.0;
        for (iy = 0; iy < FPY0_GRID; iy++)
            for (iz = 0; iz < FPZ0_GRID; iz++)

                t1 += rho[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];


        zvec[ix + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = FNX_GRID;
    global_sums(zvec, &ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < FNX_GRID; ix++)
        {
            t1 = ix * ct.hxxgrid;
            printf(" %d %f %s\n", ix, zvec[ix] / FNY_GRID / FNZ_GRID, ab);
        }
        fflush(NULL);
    }

    my_free(zvec);
}                               /* end get_avgd */

/******/
