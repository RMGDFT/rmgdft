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
#include "md.h"

#define min(a,b) (((a)>(b)) ? (b) : (a))
#define DELTA_V_MAX 1.0

void update_pot(double *, double *, double *, double *, double *, double *,
                double *, double *, int *, STATE * states);
void pulay_rho (int step0, int N, double *xm, double *fm, int NsavedSteps,
                int Nrefresh, double scale, int preconditioning);
static double t[2];
extern int it_scf;
double tem1;

void scf(STATE * states, STATE * states1, double *vxc, double *vh,
         double *vnuc, double *rho, double *rhoc, double *rhocore,
         REAL * vxc_old, REAL * vh_old, int *CONVERGENCE)
{
    double time1, time2, time3, time4;
    int numst = ct.num_states;
    int ispin, kpt, kpt1;
    int st1, idx, ione = 1;
    double tem;
    int flag;
    int steps;

    time1 = my_crtc();


    /* decide if update orbital localization centers  */
    ct.move_centers_at_this_step = 0;
    if ((ct.scf_steps % ct.movingSteps == 0) && if_update_centers(states)
        && ct.scf_steps > 1 && ct.movingCenter)
    {
        ct.move_centers_at_this_step = 1;
        update_orbit_centers(states);
        get_all_kbpsi(states, states);
        duplicate_states_info(states, states1);
        duplicate_states_info(states, states_tem);
    }


    /* Update the wavefunctions */
    time3 = my_crtc();


    if(ct.scf_steps < ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
        mg_eig(states, states1, vxc, vh, vnuc, rho, rhoc, vxc_old, vh_old);
    }
    else
    {
        steps = ct.scf_steps - ct.freeze_orbital_step;
        if(ct.charge_pulay_order ==1 )  ct.charge_pulay_order++;
    }
    time4 = my_crtc();
    rmg_timings(MG_TIME, time4 - time3);



    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        for (kpt = pct.kstart; kpt < pct.kend; kpt++)
        {
            kpt1 = kpt + ispin * ct.num_kpts;
            flag = 0;
            time3 = my_crtc();
            matrix_and_diag(ct.kp[kpt1].kstate, states1, vtot_c, flag);
            time4 = my_crtc();
            rmg_timings(MATDIAG_TIME, time4 - time3);
        }
    }

    /* Generate new density */
    ct.efermi = fill(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    scopy(&FP0_BASIS, rho, &ione, rho_old, &ione);

    get_new_rho(states, rho);

    tem1 = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
        tem1 += rho_old[idx] * rho_old[idx];
    }

    tem1 = sqrt(real_sum_all (tem1) ) /(double) FP0_BASIS;
    pulay_rho (steps, FP0_BASIS, rho, rho_old, ct.charge_pulay_order, ct.charge_pulay_refresh, ct.mix, 0); 


    /* Update potential */
    time3 = my_crtc();
    update_pot(vxc, vh, vxc_old, vh_old, vnuc, rho, rhoc, rhocore, CONVERGENCE, states);
    time4 = my_crtc();
    rmg_timings(UPDATEPOT_TIME, time4 - time3);


    get_te(rho, rhoc, rhocore, vh, vxc, states);


    time2 = my_crtc();
    rmg_timings(SCF_TIME, time2 - time1);

}                               /* end scf */



/*
   Function to update potentials vh and vxc:

   The new potentials are computed as a linear combination 
   of the old ones (input "vh" and "vxc") and the ones 
   corresponding to the input "rho".
 */
void update_pot(double *vxc, double *vh, REAL * vxc_old, REAL * vh_old,
        double *vnuc, double *rho, double *rhoc, double *rhocore,
        int *CONVERGENCE, STATE * states)
{
    int n = FP0_BASIS, idx, ione = 1;

    /* save old vtot, vxc, vh */
    scopy(&n, vxc, &ione, vxc_old, &ione);
    scopy(&n, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx];

    /* Generate exchange-correlation potential */
    get_vxc(rho, rhocore, vxc);

    pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);

    /* Generate hartree potential */
    get_vh(rho, rhoc, vh, 15, ct.poi_parm.levels);


    /* evaluate correction vh+vxc */
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] - vtot[idx];

    /* evaluate SC correction */
    t[0] = t[1] = 0.;

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        t[0] += rho[idx] * vtot[idx];
        t[1] += vtot[idx] * vtot[idx];
    }
    idx = 2;
    global_sums(t, &idx, pct.grid_comm);
    t[0] *= ct.vel_f;
    t[0] /= (double) ct.num_ions;
    t[1] = sqrt(t[1] / ((double) (ct.vh_nbasis)));

    if (pct.gridpe == 0)
        printf(" SCF CHECKS: RMS[dv] = %15.10e RMS[drho] = %15.10e \n", t[1], tem1);

    fflush(NULL);
    my_barrier();

    if (ct.scf_steps < 4 && ct.runflag == 0)
    {
        for (idx = 0; idx < FP0_BASIS; idx++)
        {
            vxc[idx] = vxc_old[idx];
            vh[idx] = vh_old[idx];
        }
    }

    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

   get_vtot_psi(vtot_c, vtot, FG_NX);
    
    if (t[1] < ct.thr_rms)
        *CONVERGENCE = TRUE;
}

