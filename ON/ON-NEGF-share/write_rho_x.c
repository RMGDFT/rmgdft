/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"

void write_rho_x(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;
    double *zvec;

    my_malloc_init( zvec, get_FNX_GRID(), double );
    /* Get this processors offset */
    poff = get_FPX_OFFSET();



    /* Zero out result vector */
    for (ix = 0; ix < get_FNX_GRID(); ix++)
        zvec[ix] = ZERO;


    /* Loop over this processor */
    for (ix = 0; ix < get_FPX0_GRID(); ix++)
    {
        t1 = 0.0;
        for (iy = 0; iy < get_FPY0_GRID(); iy++)
            for (iz = 0; iz < get_FPZ0_GRID(); iz++)

                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[ix + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNX_GRID();
    global_sums(zvec, &ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < get_FNX_GRID(); ix++)
        {
            t1 = ix * get_hxxgrid();
            printf(" %d %f %s\n", ix, zvec[ix] / get_FNY_GRID() / get_FNZ_GRID(), ab);
        }
        printf(" & %s\n", ab);
        fflush(NULL);
    }

    my_free(zvec);
}                               /* end get_avgd */

/******/
