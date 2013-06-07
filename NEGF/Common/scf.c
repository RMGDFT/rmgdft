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
#include "method.h"
#include "pmo.h"
#define min(a,b) (((a)>(b)) ? (b) : (a))


static double t[2];

void update_pot (double *, double *, double *, double *, double *, double *, double *, double *, double *,
                 int *);

extern int it_scf;



void scf (complex double * sigma_all, STATE * states, STATE * states_distribute, double *vxc,
          double *vh, double *vnuc, double *vext, double *rho, double *rhoc, double *rhocore,
          REAL * vxc_old, REAL * vh_old, REAL * vbias, int *CONVERGENCE)
{
    double time1, time2;
    int st1, st2, idx, idx1, ione = 1;
    int st11, st22;
    double tem;
    int i, j, k, jj, kk;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int j1, k1, jdiff, kdiff, iprobe, idx_C;
    int idx2, FPYZ0_GRID;


    time1 = my_crtc ();
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vh[idx] + vxc[idx] -vh_old[idx] - vxc_old[idx];



    get_vtot_psi(vtot_c, vtot, FG_NX);

    get_ddd_update (vtot);


    idx = ct.num_states - lcr[2].num_states;
    idx1 = ct.num_states - lcr[2].num_states / 2;

    /* get lcr[0].H00 part */
    get_Hij_update (states, states_distribute, vtot_c, work_matrix);


    whole_to_tri_update (lcr[0].Htri, work_matrix, ct.num_blocks, ct.block_dim);


/* ========= interaction between L3-L4 is zero ========== */

      zero_lead_image(lcr[0].Htri);  





    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    if (ct.runflag == 111 && cei.num_probe == 2)
    {
        for(j=0; j < pmo.mxllda_lead[1]; j++)
        {
            for(k =0; k < pmo.mxlocc_lead[1]; k++)
            {

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb;

                st11 = jj + lcr[1].num_states;
                st22 = kk + lcr[1].num_states;
                lcr[1].H00[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];
                lcr[2].H00[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];

                st11 = jj + lcr[1].num_states + lcr[0].num_states;
                lcr[1].H01[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];
                lcr[1].HCL[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];

                st11 = jj + lcr[1].num_states;
                st22 = kk + lcr[1].num_states + lcr[0].num_states;
                lcr[2].H01[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];
                lcr[2].HCL[j + k * pmo.mxllda_lead[1]]=
                    work_matrix[st11 * ct.num_states + st22];

            }
        }
    }

/* corner elements keep unchanged */
    setback_corner_matrix_H();  


    /* Generate new density */
    scopy (&FP0_BASIS, rho, &ione, rho_old, &ione);

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
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        tem += (rho[idx] - rho_old[idx]) * (rho[idx] - rho_old[idx]);
    }

    tem = real_sum_all (tem, pct.grid_comm);
    tem = sqrt (tem);

    if (pct.gridpe == 0)
        printf (" \nSCF CHECKS: <drho>/ion = %12.6e RMS[drho/GRID] = %12.6e\n",
                tem / ct.num_ions, tem / FP0_BASIS / NPES);


    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
    }

    pulay_rho (ct.steps, FP0_BASIS, rho, rho_old, cei.Npulaysave, cei.Npulayrefresh, cei.pulaymix,
            0);


    /* mix_rho(rho, rho_old, ct.mix, ct.steps, 1); */

    my_barrier ();


    update_pot (vxc, vh, vxc_old, vh_old, vnuc, vext, rho, rhoc, rhocore, CONVERGENCE);



    /*
     *    mix_vh(vh, vh_old, ct.mix, ct.steps, 1); 
     *    mix_vxc(vxc, vxc_old, ct.mix, ct.steps, 1);
     */

    time2 = my_crtc ();
    rmg_timings (SCF_TIME, time2 - time1);



}                               /* end scf */



/*
   Function to update potentials vh and vxc:

   The new potentials are computed as a linear combination 
   of the old ones (input "vh" and "vxc") and the ones 
   corresponding to the input "rho".
 */
void update_pot (double *vxc, double *vh, REAL * vxc_old, REAL * vh_old, double *vnuc, double *vext,
        double *rho, double *rhoc, double *rhocore, int *CONVERGENCE)
{

    int n = FP0_BASIS, idx, ione = 1;

    double time1, time2;

    /* allocate memory */

    time1 = my_crtc ();

    /* save old vtot, vxc, vh */
    scopy (&n, vxc, &ione, vxc_old, &ione);
    scopy (&n, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        vtot[idx] = vxc[idx] + vh[idx];
    }                           /* idx */

    /* Generate exchange-correlation potential */
    get_vxc(rho, rhocore, vxc);

    pack_vhstod (vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);

    /* Generate hartree potential */
//    get_vh (rho, rhoc, vh, vh_old, 15, ct.poi_parm.levels);
   get_vh_negf (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);
//   get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);



    /* check convergence */

    /* evaluate correction vh+vxc */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        vtot[idx] = vxc[idx] + vh[idx] - vtot[idx];
    }

    my_barrier ();

    if (ct.steps < 4 && ct.runflag == 0)
    {

        for (idx = 0; idx < FP0_BASIS; idx++)
        {
            vxc[idx] = vxc_old[idx];
            vh[idx] = vh_old[idx];
        }
    }

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vext[idx];
    }

    time2 = my_crtc ();
    rmg_timings (UPDATE_POT_TIME, time2 - time1);

}
