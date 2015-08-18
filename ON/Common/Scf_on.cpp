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


#include "prototypes_on.h"
#include "init_var.h"
#define DELTA_V_MAX 1.0

void update_pot(double *, double *, double *, double *, double *, double *, double *,
        double *, double *, int *, STATE * states);
void pulay_rho_on (int step0, int N, double *xm, double *fm, int NsavedSteps,
        int Nrefresh, double scale, int preconditioning);
static double t[2];
extern int it_scf;
double tem1;

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





    for (kpt = pct.kstart; kpt < ct.num_kpts; kpt+= pct.pe_kpoint)
    {
        flag = 0;
        RmgTimer *RT1 = new RmgTimer("2-SCF: matrix_and_diag");
        matrix_and_diag(ct.kp[kpt].kstate, states1, vtot_c, flag);
        
//        get_HS(states, states1, vtot_c, Hij_00, Bij_00);
        delete(RT1);
    my_barrier();
        RmgTimer *RTa = new RmgTimer("diag_ele");
        DiagElemental(ct.num_states, Hij_00, Bij_00);
    my_barrier();
        delete(RTa);
    }

    if(ct.spin_flag)
        get_opposite_eigvals( states );
    /* Generate new density */
    ct.efermi = fill(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag);

    if(pct.gridpe == 0) write_eigs(states);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        rmg_printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    dcopy(&nfp0, rho, &ione, rho_old, &ione);

    RmgTimer *RT2 = new RmgTimer("2-SCF: get_new_rho");
    get_new_rho(states, rho);
    delete(RT2);

    tem1 = 0.0;
    for (idx = 0; idx < nfp0; idx++)
    {
        tem = rho_old[idx];
        rho_old[idx] = -rho[idx] + rho_old[idx];
        rho[idx] = tem;
        tem1 += rho_old[idx] * rho_old[idx];
    }

    tem1 = sqrt(real_sum_all (tem1, pct.grid_comm) )/ ((double) (ct.vh_nbasis));
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
    pulay_rho_on (steps, nfp0, rho, rho_old, ct.charge_pulay_order, ct.charge_pulay_refresh, ct.mix, 0); 
    delete(RT3);

    if(ct.spin_flag) get_rho_oppo(rho, rho_oppo);

    /* Update potential */
    RmgTimer *RT4 = new RmgTimer("2-SCF: update_pot");
    update_pot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore, CONVERGENCE, states);
    delete(RT4);


    RmgTimer *RT5 = new RmgTimer("2-SCF: get_te");
    get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
    delete(RT5);

    /* Update the orbitals */



    if(ct.scf_steps < ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
        RmgTimer *RT6 = new RmgTimer("2-SCF: mg_eig");
        mg_eig(states, states1, vxc, vh, vnuc, rho, rhoc, vxc_old, vh_old);
        delete(RT6);
    }

    delete(RT);
}                               /* end scf */



/*
   Function to update potentials vh and vxc:

   The new potentials are computed as a linear combination 
   of the old ones (input "vh" and "vxc") and the ones 
   corresponding to the input "rho".
 */
void update_pot(double *vxc, double *vh, double * vxc_old, double * vh_old,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore,
        int *CONVERGENCE, STATE * states)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int idx, ione = 1;
    double tem;




    /* save old vtot, vxc, vh */
    dcopy(&nfp0, vxc, &ione, vxc_old, &ione);
    dcopy(&nfp0, vh, &ione, vh_old, &ione);

    for (idx = 0; idx < nfp0; idx++)
        vtot[idx] = vxc[idx] + vh[idx];

    /* Generate exchange-correlation potential */
    get_vxc(rho, rho_oppo, rhocore, vxc);

    pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.boundaryflag);

    /* Keep in memory vh*rho_new before updating vh */
    tem = ddot(&nfp0, rho, &ione, vh, &ione);
    ct.Evhold_rho = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);


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

    get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, ct.boundaryflag);



    /* Compute quantities function of rho only */
    tem = ddot(&nfp0, rho, &ione, vh, &ione);
    ct.Evh_rho = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);

    tem = ddot(&nfp0, rhoc, &ione, vh, &ione);
    ct.Evh_rhoc = 0.5 * get_vel_f() * real_sum_all(tem, pct.grid_comm);



    /* evaluate correction vh+vxc */
    for (idx = 0; idx < nfp0; idx++)
        vtot[idx] = vxc[idx] + vh[idx] - vtot[idx];


    /* evaluate SC correction */
    t[0] = t[1] = 0.;

    for (idx = 0; idx < nfp0; idx++)
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
        rmg_printf(" SCF CHECKS: RMS[dv] = %15.10e RMS[drho] = %15.10e \n", t[1], tem1);


    fflush(NULL);
    my_barrier();

    if (ct.scf_steps < 4 && ct.runflag == 0)
    {
        for (idx = 0; idx < nfp0; idx++)
        {
            vxc[idx] = vxc_old[idx];
            vh[idx] = vh_old[idx];
        }
    }

    for (idx = 0; idx < nfp0; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    get_ddd(vtot);

    get_vtot_psi(vtot_c, vtot, Rmg_G->default_FG_RATIO);

    if (t[1] < ct.thr_rms)
        *CONVERGENCE = TRUE;
}

