/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"


#if USE_DIS_MAT

#include "my_scalapack.h"

#endif


void get_Hij_update (STATE * states, STATE * states_distribute, double *vtot_c, double *Aij)
{
    int idx, st1, st2, idx1, idx2;
    int st11, st22;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    REAL tem, tem1;
    int ixx, iyy, izz;
    char msg[100];
    double time1, time2, time3, time4;
    double *psi, *mat, one = 1.0, zero = 0.0;

    int ix, iy,iz;
    time1 = my_crtc ();

    n2 = ct.num_states * ct.num_states;
    maxst = ct.num_states;


    my_malloc(psi, pct.num_local_orbit * pct.P0_BASIS+1024, double);
    my_malloc(mat, pct.num_local_orbit * pct.num_local_orbit+1024, double);
    for (st1 = 0; st1 < ct.num_states * ct.num_states; st1++)
        Aij[st1] = 0.;

    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    time3 = my_crtc ();

#if GPU_ENABLED
    cublasOperation_t transT = CUBLAS_OP_T, transN = CUBLAS_OP_N;

    cublasSetVector( pct.P0_BASIS, sizeof( double ), vtot_c, ione, ct.gpu_host_temp2, ione );
    genvpsi_gpu(ct.gpu_states, ct.gpu_host_temp2, ct.gpu_host_temp1, pct.num_local_orbit, pct.P0_BASIS);

    cublasDgemm (ct.cublas_handle, transT, transN, pct.num_local_orbit, pct.num_local_orbit, pct.P0_BASIS, 
            &one, ct.gpu_host_temp1, pct.P0_BASIS, ct.gpu_states, pct.P0_BASIS, &zero, ct.gpu_host_temp2, pct.num_local_orbit);

    cublasGetVector( pct.num_local_orbit * pct.num_local_orbit, sizeof( double ), ct.gpu_host_temp2, ione, mat, ione );

#else

    for (st1 = 0; st1 < pct.num_local_orbit; st1++)
        for(idx1 = 0; idx1 <pct.P0_BASIS; idx1++)
            psi[st1 * pct.P0_BASIS + idx1] = states_distribute[st1].psiR[idx1] * vtot_c[idx1];
    dgemm ("T", "N", &pct.num_local_orbit, &pct.num_local_orbit, &pct.P0_BASIS, &one, psi, &pct.P0_BASIS,
            states_distribute[0].psiR, &pct.P0_BASIS, &zero, mat, &pct.num_local_orbit);
#endif



    for (st1 = 0; st1 < pct.num_local_orbit; st1++)
    {
        st11 = states_distribute[st1].istate;
        for (st2 = 0; st2 < pct.num_local_orbit; st2++)
        {
            st22 = states_distribute[st2].istate;

            idx = st11 * ct.num_states + st22;
            Aij[idx] = mat[st1 * pct.num_local_orbit + st2];

        }
    }                           /* end for st1 = .. */


    my_barrier();
    time4 = my_crtc ();
    rmg_timings (H_psi_TIME, (time4 - time3));


    time3 = my_crtc ();

    get_Hvnlij (Aij);

    time4 = my_crtc ();
    rmg_timings (get_Hnl_TIME, (time4 - time3));


    /* symmetrize the Aij */
    for (st1 = 0; st1 < ct.num_states - 1; st1++)
        for (st2 = st1 + 1; st2 < ct.num_states; st2++)
        {
            idx1 = st1 + st2 * ct.num_states;
            idx2 = st2 + st1 * ct.num_states;
            Aij[idx1] = 0.5 * (Aij[idx1] + Aij[idx2]);
            Aij[idx2] = Aij[idx1];
        }

    global_sums (Aij, &n2, pct.grid_comm);     /* sum up Aij contributions */

    dscal (&n2, &ct.vel, Aij, &ione);

    if (pct.gridpe == 0)
    {
        printf (" matrix Hij\n");
        print_matrix (work_matrix, 5, maxst);
        print_sum (n2, work_matrix, "work_matrix before leaving get_Hij ddd ");
    }

    time2 = my_crtc ();

    rmg_timings (GET_Hij_TIME, (time2 - time1));

    my_free(mat);
    my_free(psi);

}
