/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"

#include "my_scalapack.h"


void get_new_rho_soft (STATE * states, double *rho)
{
    int idx, ione = 1;
    double t2;
    register double tcharge;

    /* for parallel libraries */

    double *psi1, *psi2, scale;
    int i, st1, st2, proc1, proc2, st11;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    double *rho_temp;
    char filename[MAX_PATH];

    int ix, iy;
    double tem;

    state_per_proc = ct.state_per_proc + 2;
    /*    if (pct.gridpe == 0)
          printf (" Compute new density\n");*/

    my_malloc_init( rho_global, get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(), double );
    my_malloc_init( rho_temp, get_P0_BASIS(), double );



    tri_to_row (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);

//    for(i = 0; i < pct.num_local_orbit * pct.num_local_orbit; i++) work_matrix[i] = 0.0;
//    work_matrix[100 * (pct.num_local_orbit +1)] = 10000.0;

    for (idx = 0; idx < get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(); idx++)
    {
        rho_global[idx] = 0.;
    }

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {

        st11 = st1 - ct.state_begin;
        for (st2 = st1; st2 < ct.state_end; st2++)
        {
            if (st1 == st2)
                scale = 1.0 * work_matrix[st11 * ct.num_states + st2];
            if (st1 != st2)
                scale = 2.0 * work_matrix[st11 * ct.num_states + st2];
            psi1 = states[st1].psiR;
            psi2 = states[st2].psiR;

            if (state_overlap_or_not[st11 * ct.num_states + st2 ] == 1)
                density_orbit_X_orbit (st1, st2, scale, psi1, psi2,
                        rho_global, 0, states, orbit_overlap_region);

        }
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
                    MPI_DOUBLE, proc2, i, pct.grid_comm, &mstatus);

            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {

                st11 = st1 - ct.state_begin;
                if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    psi1 = states[st1].psiR;
                    scale = 2.0 * work_matrix[st11 * ct.num_states + st2];
                    density_orbit_X_orbit (st1, st2, scale, psi1, psi2, 
                        rho_global, 0, states, orbit_overlap_region);
                }
            }
        }

        if (num_send < num_recv)
            for (i = num_send; i < num_recv; i++)
            {
                st2 = recv_from[loop * state_per_proc + i + 2];
                size2 = states[st2].size;
                MPI_Recv (psi2, size2, MPI_DOUBLE, proc2, i, pct.grid_comm, &mstatus);
                for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                {

                    st11 = st1 - ct.state_begin;
                    if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                    {
                        psi1 = states[st1].psiR;
                        scale = 2.0 * work_matrix[st11 * ct.num_states + st2];
                        density_orbit_X_orbit (st1, st2, scale, psi1, psi2, 
                        rho_global, 0, states, orbit_overlap_region);
                    }
                }
            }

        if (num_send > num_recv)
            for (i = num_recv; i < num_send; i++)
            {
                st1 = send_to[loop * state_per_proc + i + 2];
                psi1 = states[st1].psiR;
                size1 = states[st1].size;
                MPI_Send (psi1, size1, MPI_DOUBLE, proc1, i, pct.grid_comm);
            }

        my_barrier ();

    }                           /* end of loop  */

    idx = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    global_sums (rho_global, &idx, pct.grid_comm);
    global_to_distribute (rho_global, rho_temp);



    for(ix = 0; ix < get_PX0_GRID(); ix++)
    {
        tem = 0.0;
        for(iy = 0; iy < get_PY0_GRID() * get_PZ0_GRID(); iy++) tem+= rho_temp[ix * get_PY0_GRID() * get_PZ0_GRID() + iy];
        printf ("\n %d   %f  rho_oldwww ", ix, tem);
    }


    my_free(rho_global);

    mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
            get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);

    my_free(rho_temp);

    rho_augmented(rho, work_matrix, state_begin, state_end,
            num_nonlocal_ion,
            kbpsi, max_ion_nonlocal, kbpsi_comm, ionidx_allproc);


    my_barrier ();

#if  	DEBUG
    print_sum_square (get_P0_BASIS(), rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
