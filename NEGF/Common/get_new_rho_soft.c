/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "md.h"

#if NONORTHO
#if USE_DIS_MAT
#include "my_scalapack.h"
#endif
#endif


void get_new_rho_soft (STATE * states, double *rho)
{
    int idx, ione = 1;
    REAL t2;
    register double tcharge;

#if NONORTHO
#if USE_DIS_MAT
    /* for parallel libraries */
    int n2 = ct.num_states * ct.num_states;
    int mxllda;
#endif
#endif

    REAL *psi1, *psi2, scale;
    int i, st1, st2, proc1, proc2;
    REAL time1;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    REAL *rho_temp;
    char filename[MAX_PATH];


    state_per_proc = ct.state_per_proc + 2;
    time1 = my_crtc ();
/*    if (pct.gridpe == 0)
        printf (" Compute new density\n");*/

    my_malloc_init( rho_global, NX_GRID * NY_GRID * NZ_GRID, REAL );
    my_malloc_init( rho_temp, P0_BASIS, REAL );

#if  	DEBUG
    print_sum_square (P0_BASIS, rho, "rho_sum_square before get_new_rho  ");
#endif

/* #if USE_DIS_MAT && NONORTHO
 *   	mxllda = MXLLDA;
 *   	get_distributed_mat(work_matrix, mat_X);
 *   	global_sums(work_matrix, &n2);
 * #endif
 */

    tri_to_whole_p (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);


    for (idx = 0; idx < NX_GRID * NY_GRID * NZ_GRID; idx++)
    {
        rho_global[idx] = 0.;
    }

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = st1; st2 < ct.state_end; st2++)
        {
            if (st1 == st2)
                scale = 2.0 * work_matrix[st1 * ct.num_states + st2];
            if (st1 != st2)
                scale = 4.0 * work_matrix[st1 * ct.num_states + st2];
            psi1 = states[st1].psiR;
            psi2 = states[st2].psiR;

            if (state_overlap_or_not[st1 + st2 * ct.num_states] == 1)
                density_orbit_X_orbit (st1, st2, scale, psi1, psi2, rho_global, 0, states);

        }

    my_barrier ();
    psi2 = orbit_tem;
    for (loop = 0; loop < num_sendrecv_loop; loop++)
    {
        proc1 = send_to[loop * state_per_proc];
        proc2 = recv_from[loop * state_per_proc];
        num_send = send_to[loop * state_per_proc + 1];
        num_recv = recv_from[loop * state_per_proc + 1];
        num_sendrecv = min (num_send, num_recv);

        for (i = 0; i < num_sendrecv; i++)
        {
            st1 = send_to[loop * state_per_proc + i + 2];
            st2 = recv_from[loop * state_per_proc + i + 2];

            psi1 = states[st1].psiR;
            size1 = states[st1].size;
            size2 = states[st2].size;

            MPI_Sendrecv (psi1, size1, MPI_DOUBLE, proc1, i, psi2, size2,
                          MPI_DOUBLE, proc2, i, MPI_COMM_WORLD, &mstatus);

            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                {
                    psi1 = states[st1].psiR;
                    scale = 4.0 * work_matrix[st1 * ct.num_states + st2];
                    density_orbit_X_orbit (st1, st2, scale, psi1, psi2, rho_global, 0, states);
                }
        }

        if (num_send < num_recv)
            for (i = num_send; i < num_recv; i++)
            {
                st2 = recv_from[loop * state_per_proc + i + 2];
                size2 = states[st2].size;
                MPI_Recv (psi2, size2, MPI_DOUBLE, proc2, i, MPI_COMM_WORLD, &mstatus);
                for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                    if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                    {
                        psi1 = states[st1].psiR;
                        scale = 4.0 * work_matrix[st1 * ct.num_states + st2];
                        density_orbit_X_orbit (st1, st2, scale, psi1, psi2, rho_global, 0, states);
                    }
            }

        if (num_send > num_recv)
            for (i = num_recv; i < num_send; i++)
            {
                st1 = send_to[loop * state_per_proc + i + 2];
                psi1 = states[st1].psiR;
                size1 = states[st1].size;
                MPI_Send (psi1, size1, MPI_DOUBLE, proc1, i, MPI_COMM_WORLD);
            }

        my_barrier ();

    }                           /* end of loop  */

    idx = NX_GRID * NY_GRID * NZ_GRID;
    global_sums (rho_global, &idx, pct.grid_comm);
    global_to_distribute (rho_global, rho_temp);

    my_free(rho_global);

    mg_prolong_MAX10 (rho, rho_temp, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
            PX0_GRID, PY0_GRID, PZ0_GRID, FG_NX, 6);

    my_free(rho_temp);

    rho_augmented (rho, work_matrix);

    my_barrier ();
    time1 = my_crtc () - time1;
    rmg_timings (GET_NEW_RHO, time1);

#if  	DEBUG
    print_sum_square (P0_BASIS, rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
