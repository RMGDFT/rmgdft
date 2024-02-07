/************************** SVN Revision Information **************************
 **    $Id: write_rho_x.c 2325 2014-05-19 20:20:45Z luw $    **
******************************************************************************/
 


#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"

void write_rho_y(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;
    double *zvec;

    my_malloc_init( zvec, get_FNY_GRID(), double );
    /* Get this processors offset */
    poff = get_FPY_OFFSET();



    /* Zero out result vector */
    for (ix = 0; ix < get_FNY_GRID(); ix++)
        zvec[ix] = ZERO;


    /* Loop over this processor */
    for (iy = 0; iy < get_FPY0_GRID(); iy++)
    {
        t1 = 0.0;
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
            for (iz = 0; iz < get_FPZ0_GRID(); iz++)

                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[iy + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNY_GRID();
    global_sums(zvec, &ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        rmg_printf("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < get_FNY_GRID(); ix++)
        {
            t1 = ix * get_hxxgrid();
            rmg_printf(" %d %f %s\n", ix, zvec[ix] / get_FNX_GRID() / get_FNZ_GRID(), ab);
        }
        rmg_printf(" & %s\n", ab);
        fflush(NULL);
    }

    my_free(zvec);
}                               /* end get_avgd */

/******/
