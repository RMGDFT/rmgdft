/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"


#include "prototypes_on.h"
#include "init_var.h"

void UpdatePot(double *vxc, double *vh, double * vxc_old, double * vh_old,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int idx, ione = 1;




    /* save old vtot, vxc, vh */
    dcopy(&nfp0, vxc, &ione, vxc_old, &ione);
    dcopy(&nfp0, vh, &ione, vh_old, &ione);

    /* Generate exchange-correlation potential */
    get_vxc(rho, rho_oppo, rhocore, vxc);

    pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.boundaryflag);

    /* Generate hartree potential */
    //    get_vh1(rho, rhoc, vh, 15, ct.poi_parm.levels);

    if(ct.spin_flag == 1)
    {
        for (idx = 0; idx < nfp0; idx++) 
            rho_tot[idx] = rho[idx] + rho_oppo[idx];
    }
    else
    {
        for (idx = 0; idx < nfp0; idx++) 
            rho_tot[idx] = rho[idx] ;
    }

    VhDriver(rho_tot, rhoc, vh);
//    get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, ct.boundaryflag);




    /* evaluate correction vh+vxc */
    for (idx = 0; idx < nfp0; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vh_corr[idx];


    get_ddd(vtot);

    get_vtot_psi(vtot_c, vtot, Rmg_G->default_FG_RATIO);

}

