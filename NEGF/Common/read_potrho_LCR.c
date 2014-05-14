/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *
 * read_potrho_LCR.c
 *
 * This routine works with the following steps:
 * 1. Patches the potentials (and rhos) of all subsystems and store them 
 * in global array(s).
 * 2. While patching the vh is adjusted to the Fermi energy of the subsystem.
 * 3. For the 3rd, 4th probes use: pot/rho_NEGF(x,y,z) = pot/rho_ON(y,x,z).
 * 4. Interpolate data in a new fine-grid used in NEGF calculation. 
 * 5. Move the data from the global-array to the distributed array.
 *
 */



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"



void read_potrho_LCR (double *vh, double *vxc, double *rho)
{
    int fhand;
    long nbytes;
    char newname[MAX_PATH + 200];
    char msg[200];



    /* Wait until everybody gets here */
    my_barrier ();




    read_potrho (vh,  1, 0);
        if(pct.gridpe ==0) printf (" vh calc is done \n"); 


/* ===================== writing pot ======================= */


    int ix, iy, iz, idx, idx2, FPYZ0, FNXY;
    double *vtot_xyplane;
    FILE *file;
    int ii, jj, kk;

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vh[idx];




    FPYZ0 = get_FPY0_GRID() * get_FPZ0_GRID();
    FNXY = get_FNX_GRID() * get_FNY_GRID();
    my_malloc_init(vtot_xyplane, FNXY, double);


    ii = get_FPX_OFFSET();
    jj = get_FPX_OFFSET();
    kk = get_FPX_OFFSET();


    for (iz = 0; iz < get_FPZ0_GRID(); iz++)
    {
        if ((iz + kk) == (get_FNZ_GRID()/6)-1)           /* for a given iz */
        {
            for (ix = 0; ix < get_FPX0_GRID(); ix++)
            {
                for (iy = 0; iy < get_FPY0_GRID(); iy++)
                {
                    idx = ix * FPYZ0 + iy * get_FPZ0_GRID() + iz;
                    idx2 = (ix + ii) * get_FNY_GRID() + (iy + jj);
                    vtot_xyplane[idx2] = vtot[idx];
                }
            }
        }
    }

    global_sums (vtot_xyplane, &FNXY, pct.grid_comm);


    if (pct.gridpe == 0)
    {          
        file = fopen ("pot_init.dat", "w");

        for (ix = 0; ix < get_FNX_GRID(); ix++)
        {
            for (iy = 0; iy < get_FNY_GRID(); iy++)
            {
                idx = iy + ix * get_FNY_GRID();

                fprintf (file, " %d %d %f \n", ix, iy, vtot_xyplane[idx]);
            }
        }

        fclose (file);
    }

    my_free (vtot_xyplane);

/* ================================================= */






    read_potrho (vxc, 0, 1);
        if(pct.gridpe ==0) printf (" vxc calc is done \n"); 
    read_potrho (rho, 0, 2);
        if(pct.gridpe ==0) printf (" rho calc is done \n"); 





    my_barrier ();

    fflush (NULL);

}                               /* end read_data */

/* ============================================================ */
