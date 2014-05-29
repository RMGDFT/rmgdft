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


void get_new_rho_local (STATE * states_distribute, double *rho)
{
    int idx, ione = 1;
    rmg_double_t t2;
    register double tcharge;

    /* for parallel libraries */

    rmg_double_t *psi1, *psi2, scale;
    int i, st1, st2, proc1, proc2, st11;
    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    MPI_Status mstatus;
    rmg_double_t *rho_temp;
    int ix, iy; 
    double tem;
    char filename[MAX_PATH];

    double *psi, one = 1.0, zero = 0.0;
    int pbasis = get_P0_BASIS();

    psi = work_memory;
    my_malloc_init( rho_temp, pbasis, rmg_double_t );




    tri_to_local (states_distribute, lcr[0].density_matrix_tri, mat_local);

    if(pct.num_local_orbit > 0)
    {
#if GPU_ENABLED
        cublasOperation_t transN = CUBLAS_OP_N, transT = CUBLAS_OP_T;

        int n2 = pct.num_local_orbit * pct.num_local_orbit;
        cublasSetVector( n2, sizeof( double ), mat_local, ione, ct.gpu_host_temp2, ione );

        cublasDgemm (ct.cublas_handle, transN, transN, pbasis, pct.num_local_orbit, pct.num_local_orbit,
                &one, ct.gpu_states, pbasis, ct.gpu_host_temp2,
pct.num_local_orbit, &zero, ct.gpu_host_temp1, pbasis);

//        cublasGetVector( get_P0_BASIS() * pct.num_local_orbit, sizeof( double ), ct.gpu_host_temp1, ione, psi, ione );
        cublasDscal (ct.cublas_handle, pbasis, &zero, ct.gpu_host_temp2, ione);
        rho_psi_times_psi(ct.gpu_host_temp1, ct.gpu_states,
ct.gpu_host_temp2, pct.num_local_orbit, pbasis);

        cublasGetVector( pbasis, sizeof( double ), ct.gpu_host_temp2, ione, rho_temp, ione );

#else


        dgemm ("N", "N", &pbasis, &pct.num_local_orbit, &pct.num_local_orbit, &one, 
                states_distribute[0].psiR, &pbasis, mat_local, &pct.num_local_orbit, 
                &zero, psi, &pbasis);

        for(idx = 0; idx < pbasis; idx++)rho_temp[idx] = 0.0;

        for(st1 = 0; st1 < pct.num_local_orbit; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += states_distribute[st1].psiR[idx] *
psi[st1 * pbasis + idx];
#endif
    }


    mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
            get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);

    my_free(rho_temp);

    tri_to_row (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);
    rho_augmented(rho, work_matrix, state_begin, state_end,
            num_nonlocal_ion,
            kbpsi, max_ion_nonlocal, kbpsi_comm, ionidx_allproc);


    my_barrier ();

#if  	DEBUG
    print_sum_square (pbasis, rho, "rho_sum_sqare in the end of get_new_rho  ");
#endif

}
