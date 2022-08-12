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
#include "Functional.h"
#include "RmgParallelFft.h"



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


//    get_vxc(rho, rho_oppo, rhocore, vxc);

    pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.boundaryflag);

    /* Generate hartree potential */

    /* Generate exchange-correlation potential */
//double *rho_temp=new double[4*nfp0]();
//Smooth(rho_tot, rho_temp, FPX0_GRID, FPY0_GRID, FPZ0_GRID, 40.0);
//Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
//F->v_xc(rho_temp, rhocore, ct.XC, vtxc, vxc, ct.spin_flag );
//delete F;
//delete [] rho_temp;
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    F->v_xc(rho, rhocore, ct.XC, ct.vtxc, vxc, ct.nspin );
    if(ct.filter_dpot && (Rmg_G->default_FG_RATIO > 1))
    {
        for(int is = 0; is < ct.nspin; is++)
            FftFilter(&vxc[is*nfp0], *fine_pwaves, *coarse_pwaves, LOW_PASS);
    }

    delete F;


    double rms_target = ct.rms/ct.hartree_rms_ratio;
    // And new hartree potential
    VhDriver(rho, rhoc, vh, ct.vh_ext, rms_target);


    /* evaluate correction vh+vxc */
    for (idx = 0; idx < nfp0; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vh_corr[idx];

    // Transfer vtot from the fine grid to the wavefunction grid
    GetVtotPsi (vtot_c, vtot, Rmg_G->default_FG_RATIO);

    FftFilter(vtot, *fine_pwaves, *coarse_pwaves, LOW_PASS);
    get_ddd(vtot, vxc, true);


}


