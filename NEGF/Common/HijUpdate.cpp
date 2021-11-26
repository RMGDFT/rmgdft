#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "Scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

#include "pmo.h"
#include "LocalObject.h"
#include "GpuAlloc.h"



void HijUpdate (double *vtot_c)
{

    RmgTimer *RT = new RmgTimer("2-SCF: HijUpdate");

    int pbasis = Rmg_G->get_P0_BASIS(1);
    double alpha = 1.0;
    int ione = 1;

    int max_block_size = *std::max_element(ct.block_dim, ct.block_dim + ct.num_blocks);
    double *H_tem, *H_local, *H_dist;
    size_t size = max_block_size * max_block_size * sizeof(double);
    H_tem = (double *)RmgMallocHost(size);
    H_dist = (double *)RmgMallocHost(size);

    size = LocalOrbital->num_thispe * LocalOrbital->num_thispe * sizeof(double);
    H_local = (double *)RmgMallocHost(size);

    for(int st1 = 0; st1 < LocalOrbital->num_thispe; st1++)
    {
        double *a_phi = &LocalOrbital->storage_cpu[st1 * pbasis];
        double *h_phi = &H_LocalOrbital->storage_cpu[st1 * pbasis];

        for (int idx = 0; idx < pbasis; idx++)
        {
            h_phi[idx] = a_phi[idx] * vtot_c[idx];
        }
    }

    LO_x_LO(*LocalOrbital, *H_LocalOrbital, H_local, *Rmg_G);

    for(int ib = 0; ib < ct.num_blocks; ib++)
    {

        int st0 = pmo.orb_index[ib];
        int st1 = pmo.orb_index[ib+1];
        mat_local_to_glob(H_local, H_tem, *LocalOrbital, *LocalOrbital, st0, st1, st0, st1, 0);
        GetHvnlij_proj(H_tem, NULL, Kbpsi_mat_blocks[ib], Kbpsi_mat_blocks[ib],
                ct.block_dim[ib], ct.block_dim[ib], LocalProj->num_tot, false);
        int idx = (st1-st0) * (st1-st0);
        MPI_Allreduce(MPI_IN_PLACE, H_tem, idx, MPI_DOUBLE, MPI_SUM, LocalOrbital->comm);

        
        

        int *desca = &pmo.desc_cond[(ib + ib * ct.num_blocks) * DLEN];
        mat_global_to_dist(H_dist, desca, H_tem);

        int length = pmo.mxllda_cond[ib] * pmo.mxlocc_cond[ib];
        daxpy(&length, &alpha, H_dist, &ione, &lcr[0].Htri[pmo.diag_begin[ib]], &ione);

        // offdiag part
        if (ib == ct.num_blocks -1) break;
        int st2 = pmo.orb_index[ib+2];

        mat_local_to_glob(H_local, H_tem, *LocalOrbital, *LocalOrbital, st0, st1, st1, st2, 0);
        GetHvnlij_proj(H_tem, NULL, Kbpsi_mat_blocks[ib], Kbpsi_mat_blocks[ib+1],
                ct.block_dim[ib], ct.block_dim[ib+1], LocalProj->num_tot, false);

        idx = (st1-st0) * (st2-st1);
        MPI_Allreduce(MPI_IN_PLACE, H_tem, idx, MPI_DOUBLE, MPI_SUM, LocalOrbital->comm);


        desca = &pmo.desc_cond[(ib + (ib+1) * ct.num_blocks) * DLEN];
        mat_global_to_dist(H_dist, desca, H_tem);
        length = pmo.mxllda_cond[ib] * pmo.mxlocc_cond[ib+1];
        daxpy(&length, &alpha, H_dist, &ione, &lcr[0].Htri[pmo.offdiag_begin[ib]], &ione);

    }

    RmgFreeHost(H_tem);
    RmgFreeHost(H_local);
    RmgFreeHost(H_dist);

    delete(RT);

}
