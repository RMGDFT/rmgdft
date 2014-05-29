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
#include "prototypes_on.h"
#include "init_var.h"
#define min(a,b) (((a)>(b)) ? (b) : (a))
#define DELTA_V_MAX 1.0

void update_pot(double *, double *, double *, double *, double *, double *,
        double *, double *, int *, STATE * states);
void pulay_rho_on (int step0, int N, double *xm, double *fm, int NsavedSteps,
        int Nrefresh, double scale, int preconditioning);
static double t[2];
extern int it_scf;
double tem1;

void scf(STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, double *rhocore,
        rmg_double_t * vxc_old, rmg_double_t * vh_old, int *CONVERGENCE)
{
    int numst = ct.num_states;
    int ispin, kpt, kpt1;
    int idx, ione = 1;
    double tem;
    int flag;
    int steps;
    int nfp0 = get_FP0_BASIS();


    void *RT = BeginRmgTimer("2-SCF");

    ct.move_centers_at_this_step = 0;
    if ((ct.scf_steps % ct.movingSteps == 0) && if_update_centers(states)
            && ct.scf_steps > 1 && ct.movingCenter)
    {
        ct.move_centers_at_this_step = 1;
        update_orbit_centers(states);
        get_all_kbpsi(states, states, ion_orbit_overlap_region_nl, projectors, kbpsi);
        duplicate_states_info(states, states1);
        duplicate_states_info(states, states_tem);
    }





    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        for (kpt = pct.kstart; kpt < pct.kend; kpt++)
        {
            kpt1 = kpt + ispin * ct.num_kpts;
            flag = 0;
            void *RT1 = BeginRmgTimer("2-SCF: matrix_and_diag");
            matrix_and_diag(ct.kp[kpt1].kstate, states1, vtot_c, flag);
            EndRmgTimer(RT1);
        }
    }

    /* Generate new density */
    ct.efermi = fill(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    scopy(&nfp0, rho, &ione, rho_old, &ione);

    void *RT2 = BeginRmgTimer("2-SCF: get_new_rho");
    get_new_rho(states, rho);
    EndRmgTimer(RT2);

    tem1 = 0.0;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
        tem1 += rho_old[idx] * rho_old[idx];
    }

    tem1 = sqrt(real_sum_all (tem1, pct.grid_comm) ) /(double) get_FP0_BASIS();
    steps = ct.scf_steps;
    void *RT3 = BeginRmgTimer("2-SCF: pulay mix");
    pulay_rho_on (steps, get_FP0_BASIS(), rho, rho_old, ct.charge_pulay_order, ct.charge_pulay_refresh, ct.mix, 0); 
    EndRmgTimer(RT3);


    /* Update potential */
    void *RT4 = BeginRmgTimer("2-SCF: update_pot");
    update_pot(vxc, vh, vxc_old, vh_old, vnuc, rho, rhoc, rhocore, CONVERGENCE, states);
    EndRmgTimer(RT4);


    void *RT5 = BeginRmgTimer("2-SCF: get_te");
    get_te(rho, rhoc, rhocore, vh, vxc, states);
    EndRmgTimer(RT5);

    /* Update the orbitals */



    if(ct.scf_steps < ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
        void *RT6 = BeginRmgTimer("2-SCF: mg_eig");
        mg_eig(states, states1, vxc, vh, vnuc, rho, rhoc, vxc_old, vh_old);
        EndRmgTimer(RT6);
    }
    else
    {
        steps = ct.scf_steps - ct.freeze_orbital_step;
        if(ct.charge_pulay_order ==1 )  ct.charge_pulay_order++;
    }

    EndRmgTimer(RT);
}                               /* end scf */



/*
   Function to update potentials vh and vxc:

   The new potentials are computed as a linear combination 
   of the old ones (input "vh" and "vxc") and the ones 
   corresponding to the input "rho".
 */
void update_pot(double *vxc, double *vh, rmg_double_t * vxc_old, rmg_double_t * vh_old,
        double *vnuc, double *rho, double *rhoc, double *rhocore,
        int *CONVERGENCE, STATE * states)
{
    int n = get_FP0_BASIS(), idx, ione = 1;
    double tem;

    /* save old vtot, vxc, vh */
    scopy(&n, vxc, &ione, vxc_old, &ione);
    scopy(&n, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vxc[idx] + vh[idx];

    /* Generate exchange-correlation potential */
    get_vxc(rho, rho, rhocore, vxc);

    pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

    /* Keep in memory vh*rho_new before updating vh */
    tem = ddot(&n, rho, &ione, vh, &ione);
    ct.Evhold_rho = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);


    /* Generate hartree potential */
    //    get_vh1(rho, rhoc, vh, 15, ct.poi_parm.levels);
    get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, ct.boundaryflag);



    /* Compute quantities function of rho only */
    tem = ddot(&n, rho, &ione, vh, &ione);
    ct.Evh_rho = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);

    tem = ddot(&n, rhoc, &ione, vh, &ione);
    ct.Evh_rhoc = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);



    /* evaluate correction vh+vxc */
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vxc[idx] + vh[idx] - vtot[idx];

    /* evaluate SC correction */
    t[0] = t[1] = 0.;

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        t[0] += rho[idx] * vtot[idx];
        t[1] += vtot[idx] * vtot[idx];
    }
    idx = 2;
    global_sums(t, &idx, pct.grid_comm);
    t[0] *= get_vel_f();
    t[0] /= (double) ct.num_ions;
    t[1] = sqrt(t[1] / ((double) (ct.vh_nbasis)));

    ct.rms = t[1];
    if (pct.gridpe == 0)
        printf(" SCF CHECKS: RMS[dv] = %15.10e RMS[drho] = %15.10e \n", t[1], tem1);

    fflush(NULL);
    my_barrier();

    if (ct.scf_steps < 4 && ct.runflag == 0)
    {
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
        {
            vxc[idx] = vxc_old[idx];
            vh[idx] = vh_old[idx];
        }
    }

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());

    if (t[1] < ct.thr_rms)
        *CONVERGENCE = TRUE;
}

