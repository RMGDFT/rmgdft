/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"

void write_rho_x(rmg_double_t * rho, char *ab)
{

    int ix, iy, iz, poff;
    rmg_double_t t1;
    rmg_double_t *zvec;
    int pyoff, pzoff;

    my_malloc_init( zvec, get_FNX_GRID(), rmg_double_t );
    /* Get this processors offset */
    poff = pct.FPX_OFFSET;



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
            t1 = ix * ct.hxxgrid;
            printf(" %d %f %s\n", ix, zvec[ix] / get_FNY_GRID() / get_FNZ_GRID(), ab);
        }
        fflush(NULL);
    }

    my_free(zvec);
}                               /* end get_avgd */

/******/
