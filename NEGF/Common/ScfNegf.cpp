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
#include "blacs.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

#include "method.h"
#include "pmo.h"
#include "Functional.h"
#include "PulayMixing.h"



static double t[2];

void update_pot (double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                 int *);

extern int it_scf;



void ScfNegf (DoubleC *sigma_all, STATE * states, double *vxc,
          double *vh, double *vnuc, double *vext, double *rho, double *rhoc, double *rhocore, double *rho_tf,
          double * vxc_old, double * vh_old, double * vbias, int *CONVERGENCE)
{
    int st1, st2, idx, idx1, ione = 1;
    int st11, st22;
    double tem;
    int *desca;
    double one = 1.0, zero = 0.0;
    int i, j, k, jj, kk;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int j1, k1, jdiff, kdiff, iprobe, idx_C;
    int idx2, FPYZ0_GRID;

    int fpbasis;
    fpbasis = get_FP0_BASIS();

    if(ct.scf_steps == 0)
    {
        if(ct.charge_mixing_type == 0) ct.charge_pulay_order = 1;
        Pulay_rho = new PulayMixing(fpbasis, ct.charge_pulay_order, ct.charge_pulay_refresh, 
                ct.mix, ct.charge_pulay_scale, pct.grid_comm); 

    }

    RmgTimer *RT = new RmgTimer("3-SCF");


    RmgTimer *RT1 = new RmgTimer("3-SCF: Hij update");
    for (idx = 0; idx < fpbasis; idx++)
        vtot[idx] = vh[idx] + vxc[idx] -vh_old[idx] - vxc_old[idx];



    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());

    get_ddd_update (vtot);


    idx = ct.num_states - lcr[2].num_states;
    idx1 = ct.num_states - lcr[2].num_states / 2;

    /* get lcr[0].H00 part */
    HijUpdate (states, vtot_c, work_matrix);



/* ========= interaction between L3-L4 is zero ========== */

      zero_lead_image(lcr[0].Htri);  





    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    if (ct.runflag == 111 && cei.num_probe == 2)
    {
        idx1 = 2*pmo.mxllda_lead[1] * pmo.mxlocc_lead[1];
        idx2 = 3*pmo.mxllda_lead[1] * pmo.mxlocc_lead[1];
 
        for(j=0; j < pmo.mxllda_lead[1]; j++)
        {
            for(k =0; k < pmo.mxlocc_lead[1]; k++)
            {

                lcr[1].H00[j + k * pmo.mxllda_lead[1]]=
                    lcr[0].Htri[idx1 + j + k * pmo.mxllda_lead[1]];
                lcr[2].H00[j + k * pmo.mxllda_lead[1]]=
                    lcr[0].Htri[idx1 + j + k * pmo.mxllda_lead[1]];

                lcr[2].H01[j + k * pmo.mxllda_lead[1]]=
                    lcr[0].Htri[idx2 + j + k * pmo.mxllda_lead[1]];
                lcr[2].HCL[j + k * pmo.mxllda_lead[1]]=
                    lcr[0].Htri[idx2 + j + k * pmo.mxllda_lead[1]];


            }
        }

        int numst = lcr[1].num_states;
        desca = pmo.desc_lead;
        pdtran(&numst, &numst, &one, lcr[2].H01, &ione, &ione, desca,
                &zero, lcr[1].H01, &ione, &ione, desca);
        pdtran(&numst, &numst, &one, lcr[2].HCL, &ione, &ione, desca,
                &zero, lcr[1].HCL, &ione, &ione, desca);

    }

    /* corner elements keep unchanged */
    setback_corner_matrix_H();  
 
    delete(RT1);


    /* Generate new density */
    dcopy (&fpbasis, rho, &ione, rho_old, &ione);

    my_barrier ();
    if (ct.runflag == 111 && ct.metal == 1)
        find_fermi (sigma_all);
    if (ct.runflag == 111 && ct.metal == 0)
    {

        RmgTimer *RT2 = new RmgTimer("3-SCF: sigma_all for 3lead");

        sigma_all_energy_point (sigma_all, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);
        delete(RT2);
    }

    my_barrier ();

    RmgTimer *RT3 = new RmgTimer("3-SCF: charge_density_matrix");
    charge_density_matrix_p (sigma_all);


    my_barrier ();

    delete(RT3);

#if DEBUG |1
    write_rho_x (rho, "rhoooo_1");
    if (pct.imgpe == 0)
        rmg_printf ("\n rhoooo");
    write_rho_x (vtot, "vtot_1");
    if (pct.imgpe == 0)
        rmg_printf ("\n  vtot");
    write_rho_x (vh, "vhhh_1");
    if (pct.imgpe == 0)
        rmg_printf ("\n  vhhh");
    write_rho_x (vxc, "vxc_1");
    if (pct.imgpe == 0)
        rmg_printf ("\n  vxccch");
#endif


    //    get_new_rho_soft (states, rho);
    RmgTimer *RT4 = new RmgTimer("3-SCF: rho");
//    get_new_rho_soft (states, rho);
    tri_to_row (lcr[0].density_matrix_tri, work_matrix, ct.num_blocks, ct.block_dim);
    GetNewRho_on(states, rho, work_matrix);
//    get_new_rho_local (states_distribute, rho);
    delete(RT4);

#if DEBUG 
    write_rho_x (rho, "rhoaaa_1");
    if (pct.imgpe == 0)
        rmg_printf ("\n %rhoaaa");
#endif

    RmgTimer *RT5 = new RmgTimer("3-SCF: rho mixing");
    modify_rho (rho, rho_old); 
    /* modify_rho_y (rho, rho_old); */

    tem = 0.0;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        tem += (rho[idx] - rho_old[idx]) * (rho[idx] - rho_old[idx]);
    }

    tem = real_sum_all (tem, pct.grid_comm);
    tem = sqrt (tem);

    if (pct.imgpe == 0)
        rmg_printf (" \nSCF CHECKS: <drho>/ion = %12.6e RMS[drho/GRID] = %12.6e\n",
                tem / ct.num_ions, tem / get_FP0_BASIS() / pct.grid_npes);


    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = rho[idx] - rho_old[idx];
        rho[idx] = tem;
    }

    Pulay_rho->Mixing(rho, rho_old);


    /* mix_rho(rho, rho_old, ct.mix, ct.steps, 1); */

    my_barrier ();

    delete(RT5);

    RmgTimer *RT6 = new RmgTimer("3-SCF: pot update");
    update_pot (vxc, vh, vxc_old, vh_old, vnuc, vext, rho, rhoc, rhocore, rho_tf, CONVERGENCE);
    delete(RT6);
    delete(RT);



    /*
     *    mix_vh(vh, vh_old, ct.mix, ct.steps, 1); 
     *    mix_vxc(vxc, vxc_old, ct.mix, ct.steps, 1);
     */




}                               /* end scf */



/*
   Function to update potentials vh and vxc:

   The new potentials are computed as a linear combination 
   of the old ones (input "vh" and "vxc") and the ones 
   corresponding to the input "rho".
   */
void update_pot (double *vxc, double *vh, double * vxc_old, double * vh_old, double *vnuc, double *vext,
        double *rho, double *rhoc, double *rhocore, double *rho_tf, int *CONVERGENCE)
{

    int n = get_FP0_BASIS(), idx, ione = 1;


    /* allocate memory */


    /* save old vtot, vxc, vh */
    dcopy (&n, vxc, &ione, vxc_old, &ione);
    dcopy (&n, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {

        vtot[idx] = vxc[idx] + vh[idx];
    }                           /* idx */

    /* Generate exchange-correlation potential */
    //get_vxc(rho, rho, rhocore, vxc);
    double vtxc, etxc;
    RmgTimer *RT1 = new RmgTimer("2-Init: exchange/correlation");
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.spin_flag );
    delete F;
    delete RT1;


    pack_vhstod (vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

    if (ct.num_tfions > 0)
    {
        /*Add charge density from simplified waters to rho and get Hartree potential*/
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
        {
            rho[idx] += rho_tf[idx];
        }
    }	


    /* Generate hartree potential */
    //    get_vh (rho, rhoc, vh, vh_old, 15, ct.poi_parm.levels);
    get_vh_negf (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);
    //   get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);


    if (ct.num_tfions > 0)
    {
        /*Back out density due to simplified waters after Hartree potential is evaluated*/
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
        {
            rho[idx] -= rho_tf[idx];
        }
    }	


    /* check convergence */

    /* evaluate correction vh+vxc */
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        vtot[idx] = vxc[idx] + vh[idx] - vtot[idx];
    }

    my_barrier ();

    if (ct.scf_steps < 4 && ct.runflag == 0)
    {

        for (idx = 0; idx < get_FP0_BASIS(); idx++)
        {
            vxc[idx] = vxc_old[idx];
            vh[idx] = vh_old[idx];
        }
    }

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vext[idx];
    }


}

