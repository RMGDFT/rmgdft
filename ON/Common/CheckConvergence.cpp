/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*
   scf.c
   Performs a single self consistent step.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "main.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"
#include "rmg_sum_all.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "PulayMixing.h"
#include "LocalObject.h"
#include "Kbpsi.h"
#include "LdaU_on.h"
#define DELTA_V_MAX 1.0


void CheckConvergence(double *vxc, double *vh, double * vxc_old, double * vh_old, double *rho, double *rho_pre, int *CONVERGENCE)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int idx;
    double drho_max, dvh_max, dvxc_max;


    //  maximum value in rho_old = rho-rho_prev_step
    for (idx = 0; idx < nfp0; idx++)
        rho_old[idx] = fabs(rho[idx] - rho_pre[idx]);
    drho_max = *std::max_element(rho_old, rho_old+nfp0);
    //    idx = idamax(&nfp0, rho_old, &ione);
    //    drho_max = fabs(rho_old[idx]);

    for (idx = 0; idx < nfp0; idx++)
        rho_old[idx] = fabs(vh[idx] - vh_old[idx]);
    dvh_max = *std::max_element(rho_old, rho_old+nfp0);
    //    idx = idamax(&nfp0, rho_old, &ione);
    //    dvh_max = fabs(rho_old[idx]);

    dvxc_max = 0.0;

    for (idx = 0; idx < nfp0; idx++)
    {
        rho_old[idx] = fabs(vxc[idx] - vxc_old[idx]);
        if(rho_old[idx] > dvxc_max ) {
            dvxc_max = rho_old[idx];
        }
    }
    //dvxc_max = *std::max_element(rho_old, rho_old+nfp0);
    //    idx = idamax(&nfp0, rho_old, &ione);
    //    dvxc_max = fabs(rho_old[idx]);


    drho_max = rmg::max_all<double>(drho_max, pct.grid_comm); 
    dvh_max = rmg::max_all<double>(dvh_max, pct.grid_comm); 
    dvxc_max = rmg::max_all<double>(dvxc_max, pct.grid_comm); 


    ct.rms = fabs(dvh_max) + fabs(dvxc_max);

    *CONVERGENCE = FALSE;
    if (pct.gridpe == 0)
    {
        rmg::printlog(" \nSCF CHECKS: RMS[maxdvh ] = %10.5E", dvh_max);
        rmg::printlog(" [maxdrho] = %10.5E", drho_max);
        rmg::printlog(" [maxdvxc] = %10.5E\n", dvxc_max);


        fflush(NULL);

    }
    if (ct.rms < ct.thr_rms)
        *CONVERGENCE = TRUE;

}

