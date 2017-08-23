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
#include "RmgSumAll.h"


#include "prototypes_on.h"
#include "init_var.h"
#define DELTA_V_MAX 1.0

void update_pot(double *, double *, double *, double *, double *, double *, double *,
        double *, double *, int *, STATE * states);
void pulay_rho_on (int step0, int N, double *xm, double *fm, int NsavedSteps,
        int Nrefresh, double scale, int preconditioning);
extern int it_scf;

void Scf_on(STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore,
        double * vxc_old, double * vh_old, int *CONVERGENCE)
{
    int numst = ct.num_states;
    int  kpt;
    int idx, ione = 1;
    double tem;
    int flag;
    int steps;
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int px0_grid = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    double *rho_pre;

    rho_pre = new double[nfp0];
    RmgTimer *RT = new RmgTimer("2-SCF");

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


    flag = 0;
    my_barrier();



    RmgTimer *RT0 = new RmgTimer("2-SCF: orbital_comm");
    orbital_comm(states);
    delete(RT0);

    RmgTimer *RTk = new RmgTimer("2-SCF: kbpsi");
    KbpsiUpdate(states);
    delete(RTk);

    RmgTimer *RT1 = new RmgTimer("2-SCF: get_HS");
    GetHS(states, states1, vtot_c, Hij_00, Bij_00);
    delete(RT1);
#if ELEMENTAL_LIBS
    RmgTimer *RTa = new RmgTimer("2-SCF: DiagElemental");
    DiagElemental(states, ct.num_states, Hij_00, Bij_00, work_matrix_row, theta);
    delete(RTa);
#else
    RmgTimer *RTb = new RmgTimer("2-SCF: DiagScalapack");
    DiagScalapack(states, ct.num_states, Hij_00, Bij_00, work_matrix_row, theta);
    delete(RTb);
#endif


    if(ct.spin_flag)
    {
        printf("\n need to split Diag* to take care of fill for spin_polarized");
        fflush(NULL);
        exit(0);
        get_opposite_eigvals( states );
    }
    /* Generate new density */
    //    ct.efermi = fill(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag);

    if(pct.gridpe == 0) write_eigs(states);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        rmg_printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    dcopy(&nfp0, rho, &ione, rho_old, &ione);
    dcopy(&nfp0, rho, &ione, rho_pre, &ione);

    RmgTimer *RT2 = new RmgTimer("2-SCF: get_new_rho");
    GetNewRho_on(states, rho, work_matrix_row);
    //BroydenPotential(rho_old, rho, rhoc, vh_old, vh, ct.charge_broyden_order, false);

    int iii = get_FP0_BASIS();

    double tcharge = 0.0;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        tcharge += rho[idx];
    ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
    ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


    ct.tcharge *= get_vel_f();

    double t2 = ct.nel / ct.tcharge;
    dscal(&iii, &t2, &rho[0], &ione);


    if(fabs(t2 -1.0) > 1.0e-6 && pct.gridpe == 0)
        printf("\n Warning: total charge Normalization constant = %e  \n", t2);


    delete(RT2);

    for (idx = 0; idx < nfp0; idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
    }

    RmgTimer *RT3 = new RmgTimer("2-SCF: pulay mix");
    if(ct.scf_steps <ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
    }
    else
    {
        if(ct.charge_pulay_order ==1 )  ct.charge_pulay_order++;
        steps = ct.scf_steps - ct.freeze_orbital_step;
    }
    get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
    pulay_rho_on (steps, nfp0, rho, rho_old, ct.charge_pulay_order, ct.charge_pulay_refresh, ct.mix, 0); 
//for(int i=0;i<nfp0;i++)rho[i] = (1.0-ct.mix)*rho_pre[i] + ct.mix*rho[i];
    delete(RT3);

    if(ct.spin_flag) get_rho_oppo(rho, rho_oppo);

    /* Update potential */
    RmgTimer *RT4 = new RmgTimer("2-SCF: update_pot");
    UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);
    delete(RT4);

    CheckConvergence(vxc, vh, vxc_old, vh_old, rho, rho_pre, CONVERGENCE);

    RmgTimer *RT5 = new RmgTimer("2-SCF: get_te");
    delete(RT5);

    /* Update the orbitals */



    if(ct.scf_steps < ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
        RmgTimer *RT6 = new RmgTimer("2-SCF: OrbitalOptimize");
        OrbitalOptimize(states, states1, vxc, vh, vnuc, rho, rhoc, vxc_old, vh_old);
        delete(RT6);
    }

    delete(RT);

    delete [] rho_pre;
}                               /* end scf */



void CheckConvergence(double *vxc, double *vh, double * vxc_old, double * vh_old, double *rho, double *rho_pre, int *CONVERGENCE)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int idx, ione = 1;
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
    int idx_max = -100;

    for (idx = 0; idx < nfp0; idx++)
    {
        rho_old[idx] = fabs(vxc[idx] - vxc_old[idx]);
        if(rho_old[idx] > dvxc_max ) {
            dvxc_max = rho_old[idx];
            idx_max = idx;
        }
    }
    //dvxc_max = *std::max_element(rho_old, rho_old+nfp0);
//    idx = idamax(&nfp0, rho_old, &ione);
//    dvxc_max = fabs(rho_old[idx]);

    
    drho_max = RmgMaxAll<double>(drho_max, pct.grid_comm); 
    dvh_max = RmgMaxAll<double>(dvh_max, pct.grid_comm); 
    dvxc_max = RmgMaxAll<double>(dvxc_max, pct.grid_comm); 


    ct.rms = fabs(dvh_max) + fabs(dvxc_max);
    
    *CONVERGENCE = FALSE;
    if (pct.gridpe == 0)
    {
        rmg_printf(" \nSCF CHECKS: RMS[maxdvh ] = %10.5E", dvh_max);
        rmg_printf(" [maxdrho] = %10.5E", drho_max);
        rmg_printf(" [maxdvxc] = %10.5E\n", dvxc_max);


        fflush(NULL);

    }
    if (ct.rms < ct.thr_rms)
        *CONVERGENCE = TRUE;

}
