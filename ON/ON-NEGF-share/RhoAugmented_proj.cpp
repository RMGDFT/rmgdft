/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"

#include "prototypes_on.h"
#include "LocalObject.h"
#include "GpuAlloc.h"
#include "RmgGemm.h"
#include "Kbpsi.h"


// rho_aug = Qnm(r) <beta_n|phi_i> rho_mat_ij <phi_j|beta_m>

void RhoAugmented_proj(double * rho, double *rho_matrix_local)

{

    int num_orb = LocalOrbital->num_thispe;
    int num_prj = LocalProj->num_thispe;

    if(num_orb < 1 || num_prj < 1) return;

    double *rho_kbpsi = (double *)RmgMallocHost(num_orb * num_prj * sizeof(double));
    double *Qnm_coeff = (double *)RmgMallocHost(num_prj * num_prj * sizeof(double));

    RmgTimer *RT5 = new RmgTimer("3-get_new_rho: Qnm_coeff)");
    double one(1.0), zero(0.0);
    RmgGemm("N", "C", num_orb, num_prj, num_orb, one, rho_matrix_local, num_orb,
                Kbpsi_mat_local, num_prj, zero, rho_kbpsi, num_orb);

    // don't need the whole Qnm_coeff matrix, only need the diagonal blocks for ions, each block is nh * nh
    // if this is slow, we can split and only calculate the necessary blocks
    RmgGemm("N", "N", num_prj, num_prj, num_orb, one, Kbpsi_mat_local, num_prj,
                rho_kbpsi, num_orb, zero, Qnm_coeff, num_prj);
    delete RT5;

    int idx;
    int *ivec,  idx1, idx2;
    int nh, icount, ncount, i, j, ion;
    ION *iptr;
    SPECIES *sp;


    RmgTimer *RT6 = new RmgTimer("3-get_new_rho: augmented_Q(r)");
    int proj_index = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];
        nh = sp->num_projectors;
        int proj_local_index = LocalProj->index_global_to_proj[proj_index];
        proj_index += nh;
        if (proj_local_index == -1) 
        {
            //  this ion has no overlap with this processor
            continue;
        }

        ivec = Atoms[ion].Qindex.data();
        ncount = Atoms[ion].Qindex.size();

        if (Atoms[ion].Qindex.size())
        {

            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    idx1 = (proj_local_index + i) * num_prj + proj_local_index + j;
                    idx2 = (proj_local_index + j) * num_prj + proj_local_index + i;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        double Qr = GetAugcharge(i, j, icount, ct.cg_coeff.data(), iptr);
                        if (i != j)
                            rho[ivec[icount]] +=
                                Qr * (Qnm_coeff[idx1] + Qnm_coeff[idx2]);
                        else
                            rho[ivec[icount]] += Qr * Qnm_coeff[idx1];
                    }           /*end for icount */
                    idx++;
                }               /*end for j */
            }                   /*end for i */

        }                       /*end if */

    }                           /*end for ion */

    delete(RT6);
    /* release our memory */

    RmgFreeHost(rho_kbpsi);
    RmgFreeHost(Qnm_coeff);
}

