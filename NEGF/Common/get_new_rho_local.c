/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"

#include "my_scalapack.h"


void get_new_rho_local (STATE * states_distribute, double *rho)
{
    int idx, ione = 1;
    REAL t2;
    register double tcharge;

    /* for parallel libraries */

    REAL *psi1, *psi2, scale;
    int i, st1, st2, proc1, proc2, st11;
    REAL time1;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    REAL *rho_temp;
    int ix, iy; 
    double tem;
    char filename[MAX_PATH];

    double *psi, one = 1.0, zero = 0.0;
    time1 = my_crtc ();

    my_malloc(psi, pct.num_local_orbit * pct.P0_BASIS+1024, double);
    my_malloc_init( rho_temp, P0_BASIS, REAL );




    tri_to_local (states_distribute, lcr[0].density_matrix_tri, work_matrix);

    dgemm ("N", "N", &pct.P0_BASIS, &pct.num_local_orbit, &pct.num_local_orbit, &one, 
            states_distribute[0].psiR, &pct.P0_BASIS, work_matrix, &pct.num_local_orbit, 
            &zero, psi, &pct.P0_BASIS);

    for(idx = 0; idx < pct.P0_BASIS; idx++)rho_temp[idx] = 0.0;

    for(st1 = 0; st1 < pct.num_local_orbit; st1++)
        for(idx = 0; idx < pct.P0_BASIS; idx++)
            rho_temp[idx] += states_distribute[st1].psiR[idx] * psi[st1 * pct.P0_BASIS + idx];


    for(ix = 0; ix < pct.PX0_GRID; ix++)
    {
        tem = 0.0;
        double tem1 = 0.0;
        double tem2 = 0.0;
        for(iy = 0; iy < pct.PY0_GRID * pct.PZ0_GRID; iy++) tem+= rho_temp[ix * pct.PY0_GRID * pct.PZ0_GRID + iy];
        for(iy = 0; iy < pct.PY0_GRID * pct.PZ0_GRID; iy++) tem1+= states_distribute[100].psiR[ix * pct.PY0_GRID * pct.PZ0_GRID + iy]*states_distribute[100].psiR[ix * pct.PY0_GRID * pct.PZ0_GRID + iy];
        for(iy = 0; iy < pct.PY0_GRID * pct.PZ0_GRID; iy++) tem2+= psi[100 * pct.P0_BASIS+ ix * pct.PY0_GRID * pct.PZ0_GRID + iy]*psi[100 * pct.P0_BASIS+ ix * pct.PY0_GRID * pct.PZ0_GRID + iy];
        printf ("\n %d   %f  %f  %frho_newwww ", ix, tem, tem1, tem2);
    }


    mg_prolong_MAX10 (rho, rho_temp, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
            PX0_GRID, PY0_GRID, PZ0_GRID, FG_NX, 6);

    my_free(rho_temp);
    my_free(psi);

    tri_to_row (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);
    rho_augmented (rho, work_matrix);

    my_barrier ();
    time1 = my_crtc () - time1;
    rmg_timings (GET_NEW_RHO, time1);

#if  	DEBUG
    print_sum_square (P0_BASIS, rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
