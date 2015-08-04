/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"



#include "my_scalapack.h"


void get_Hij_update (STATE * states, STATE * states_distribute, double *vtot_c, double *Aij)
{
    int idx, st1, st2, idx1, idx2;
    int st11, st22;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    rmg_double_t tem, tem1;
    int ixx, iyy, izz;
    char msg[100];
    double *psi, one = 1.0, zero = 0.0;

    int ix, iy,iz;

    maxst = ct.num_states;
    int pbasis = get_P0_BASIS();


    my_malloc(psi, pct.num_local_orbit * get_P0_BASIS()+1024, double);
    for (st1 = 0; st1 < ct.num_states * (ct.state_end-ct.state_begin); st1++)
        Hij_00[st1] = 0.;

    get_Hvnlij (Hij_00, Bij_00);
    double vel = get_vel();

    n2 = ct.num_states * (ct.state_end-ct.state_begin);
    dscal (&n2, &vel, Hij_00, &ione);
    row_to_tri_p (lcr[0].Htri, Hij_00, ct.num_blocks, ct.block_dim);

    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */


    if(pct.num_local_orbit >0)
    {
//#if GPU_ENABLED
//    cublasOperation_t transT = CUBLAS_OP_T, transN = CUBLAS_OP_N;
//
//    cublasSetVector( get_P0_BASIS(), sizeof( double ), vtot_c, ione, ct.gpu_host_temp2, ione );
//    genvpsi_gpu(ct.gpu_states, ct.gpu_host_temp2, ct.gpu_host_temp1, pct.num_local_orbit, get_P0_BASIS());
//
//    cublasDgemm (ct.cublas_handle, transT, transN, pct.num_local_orbit, pct.num_local_orbit, get_P0_BASIS(), 
 //           &one, ct.gpu_host_temp1, get_P0_BASIS(), ct.gpu_states, get_P0_BASIS(), &zero, ct.gpu_host_temp2, pct.num_local_orbit);

//    cublasGetVector( pct.num_local_orbit * pct.num_local_orbit, sizeof( double ), ct.gpu_host_temp2, ione, mat_local, ione );
//
//#else

    for (st1 = 0; st1 < pct.num_local_orbit; st1++)
        for(idx1 = 0; idx1 <get_P0_BASIS(); idx1++)
            psi[st1 * get_P0_BASIS() + idx1] = states_distribute[st1].psiR[idx1] * vtot_c[idx1];
    dgemm ("T", "N", &pct.num_local_orbit, &pct.num_local_orbit,
&pbasis, &one, psi, &pbasis,
            states_distribute[0].psiR, &pbasis, &zero, mat_local, &pct.num_local_orbit);
//#endif

    }


    my_barrier();

    n2 = pct.num_local_orbit * pct.num_local_orbit;
    dscal (&n2, &vel, mat_local, &ione);

    local_to_tri(states_distribute, lcr[0].Htri, mat_local);




    my_free(psi);

}
