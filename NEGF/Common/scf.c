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
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "method.h"
#include "pmo.h"
#define min(a,b) (((a)>(b)) ? (b) : (a))


static double t[2];

void update_pot (double *, double *, double *, double *, double *, double *, double *, double *, double *,
                 int *);

extern int it_scf;



void scf (complex double * sigma_all, STATE * states, STATE * states_distribute, double *vxc,
          double *vh, double *vnuc, double *vext, double *rho, double *rhoc, double *rhocore,
          rmg_double_t * vxc_old, rmg_double_t * vh_old, rmg_double_t * vbias, int *CONVERGENCE)
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

    for (idx = 0; idx < fpbasis; idx++)
        vtot[idx] = vh[idx] + vxc[idx] -vh_old[idx] - vxc_old[idx];



    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());

    get_ddd_update (vtot);


    idx = ct.num_states - lcr[2].num_states;
    idx1 = ct.num_states - lcr[2].num_states / 2;

    /* get lcr[0].H00 part */
    get_Hij_update (states, states_distribute, vtot_c, work_matrix);



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
        PDTRAN(&numst, &numst, &one, lcr[2].H01, &ione, &ione, desca,
                &zero, lcr[1].H01, &ione, &ione, desca);
        PDTRAN(&numst, &numst, &one, lcr[2].HCL, &ione, &ione, desca,
                &zero, lcr[1].HCL, &ione, &ione, desca);

    }

    /* corner elements keep unchanged */
    setback_corner_matrix_H();  


    /* Generate new density */
    scopy (&fpbasis, rho, &ione, rho_old, &ione);

    my_barrier ();
    if (ct.runflag == 111 && ct.metal == 1)
        find_fermi (sigma_all);
    if (ct.runflag == 111 && ct.metal == 0)
        sigma_all_energy_point (sigma_all);


    my_barrier ();

    charge_density_matrix_p (sigma_all);


    my_barrier ();


#if DEBUG |1
    write_rho_x (rho, "rhoooo_1");
    if (pct.gridpe == 0)
        printf ("\n rhoooo");
    write_rho_x (vtot, "vtot_1");
    if (pct.gridpe == 0)
        printf ("\n  vtot");
    write_rho_x (vh, "vhhh_1");
    if (pct.gridpe == 0)
        printf ("\n  vhhh");
    write_rho_x (vxc, "vxc_1");
    if (pct.gridpe == 0)
        printf ("\n  vxccch");
#endif


    //    get_new_rho_soft (states, rho);
    get_new_rho_local (states_distribute, rho);
#if DEBUG
    write_rho_x (rho, "rhoaaa_1");
    if (pct.gridpe == 0)
        printf ("\n %rhoaaa");
#endif

    modify_rho (rho, rho_old); 
    /* modify_rho_y (rho, rho_old); */

    tem = 0.0;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        tem += (rho[idx] - rho_old[idx]) * (rho[idx] - rho_old[idx]);
    }

    tem = real_sum_all (tem, pct.grid_comm);
    tem = sqrt (tem);

    if (pct.gridpe == 0)
        printf (" \nSCF CHECKS: <drho>/ion = %12.6e RMS[drho/GRID] = %12.6e\n",
                tem / ct.num_ions, tem / get_FP0_BASIS() / NPES);


    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
    }

    pulay_rho_on (ct.scf_steps, fpbasis, rho, rho_old, cei.Npulaysave, cei.Npulayrefresh, cei.pulaymix,
            0);


    /* mix_rho(rho, rho_old, ct.mix, ct.steps, 1); */

    my_barrier ();


    update_pot (vxc, vh, vxc_old, vh_old, vnuc, vext, rho, rhoc, rhocore, CONVERGENCE);



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
void update_pot (double *vxc, double *vh, rmg_double_t * vxc_old, rmg_double_t * vh_old, double *vnuc, double *vext,
        double *rho, double *rhoc, double *rhocore, int *CONVERGENCE)
{

    int n = get_FP0_BASIS(), idx, ione = 1;


    /* allocate memory */


    /* save old vtot, vxc, vh */
    scopy (&n, vxc, &ione, vxc_old, &ione);
    scopy (&n, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {

        vtot[idx] = vxc[idx] + vh[idx];
    }                           /* idx */

    /* Generate exchange-correlation potential */
    get_vxc(rho, rho, rhocore, vxc);

    pack_vhstod (vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

    /* Generate hartree potential */
    //    get_vh (rho, rhoc, vh, vh_old, 15, ct.poi_parm.levels);
    get_vh_negf (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);
    //   get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);



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
