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
#include "PulayMixing.h"
#include "LocalObject.h"
#include "Kbpsi.h"
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
    int idx, ione = 1;
    double tem;
    int steps;
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    double *rho_pre;

    if(ct.scf_steps == 0)
    {
        Pulay_rho = new PulayMixing(nfp0, ct.charge_pulay_order, ct.charge_pulay_refresh, 
                ct.mix, ct.charge_pulay_scale, pct.grid_comm); 
        Pulay_orbital = new PulayMixing(pct.psi_size, ct.orbital_pulay_order, ct.orbital_pulay_refresh, 
                ct.orbital_pulay_mixfirst, ct.orbital_pulay_scale, pct.grid_comm); 
        Pulay_orbital->SetPrecond(Precond);
    }
    rho_pre = new double[nfp0];
    double *trho = new double[nfp0];
    RmgTimer *RT = new RmgTimer("2-SCF");

    ct.move_centers_at_this_step = 0;
    if ((ct.scf_steps % ct.movingSteps == 0) && if_update_centers(states)
            && ct.scf_steps > 1 && ct.movingCenter)
    {
        ct.move_centers_at_this_step = 1;
        update_orbit_centers(states);
        GetAllKbpsi(states, states, ion_orbit_overlap_region_nl, projectors, kbpsi);
        duplicate_states_info(states, states1);
        duplicate_states_info(states, states_tem);
    }


    MPI_Barrier(pct.img_comm);

    RmgTimer *RT0 = new RmgTimer("2-SCF: orbital_comm");
    for (int st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        int ixx = states[st1].orbit_nx;
        int iyy = states[st1].orbit_ny;
        int izz = states[st1].orbit_nz;

        ZeroBoundary(states[st1].psiR, ixx, iyy, izz);
    }
    
    if(ct.LocalizedOrbitalLayout == LO_projection)
    {
        write_data (ct.outfile, vh, vxc, vh_old, vxc_old, rho, vh_corr, states);
        MPI_Barrier(pct.img_comm);
        LocalOrbital->ReadOrbitals(std::string(ct.outfile), Rmg_G);
    }

    OrbitalComm(states);
    delete(RT0);


    if(ct.LocalizedOrbitalLayout == LO_projection)
    {
        LO_x_LO(*LocalProj, *LocalOrbital, Kbpsi_mat, *Rmg_G);
        GetHS_dis(*LocalOrbital, *H_LocalOrbital, vtot_c, Hij, matB, Kbpsi_mat);
    }
    else
    {
        RmgTimer *RTk = new RmgTimer("2-SCF: kbpsi");
        KbpsiUpdate(states);
        delete(RTk);
        RmgTimer *RT1 = new RmgTimer("2-SCF: get_HS");
        GetHS(states, states1, vtot_c, Hij_00, Bij_00);
        MyCpdgemr2d(numst,numst, Hij_00, pct.descb, Hij, pct.desca);
        MyCpdgemr2d(numst,numst, Bij_00, pct.descb, matB, pct.desca);
        delete(RT1);
    }

    RmgTimer *RTb = new RmgTimer("2-SCF: DiagScalapack");

    DiagScalapack(states, ct.num_states, Hij, matB);

    MyCpdgemr2d(numst, numst, mat_X, pct.desca, work_matrix_row, pct.descb);
    MyCpdgemr2d(numst,numst, uu_dis, pct.desca, theta, pct.descb);
    delete(RTb);

    if(ct.spin_flag)
    {
        get_opposite_eigvals( states );
    }
    /* Generate new density */

    ct.efermi = Fill_on(states, ct.occ_width, ct.nel, ct.occ_mix, numst, ct.occ_flag, ct.mp_order);

    get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
    if(pct.gridpe == 0) write_eigs(states);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        rmg_printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    dcopy(&nfp0, rho, &ione, rho_old, &ione);
    dcopy(&nfp0, rho, &ione, rho_pre, &ione);

    RmgTimer *RT2 = new RmgTimer("2-SCF: get_new_rho");
    for (idx = 0; idx < nfp0; idx++)trho[idx] = rho[idx];

    if(ct.scf_steps >= ct.freeze_rho_steps)
    {
        if(ct.LocalizedOrbitalLayout == LO_projection)
            GetNewRho_dis(*LocalOrbital, *H_LocalOrbital, rho, work_matrix_row);
        else
            GetNewRho_on(states, rho, work_matrix_row);
    }
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

    if(ct.spin_flag)
        get_rho_oppo(rho, rho_oppo);

    //    if(fabs(t2 -1.0) > 1.0e-6 && pct.gridpe == 0)
    printf("\n Warning: total charge Normalization constant = %15.12e  \n", t2);
    delete(RT2);

    RmgTimer *RT3 = new RmgTimer("2-SCF: pulay mix");

    if(ct.charge_mixing_type == 0)
    {
        if(ct.spin_flag)
            get_rho_oppo(rho, rho_oppo);
        //     get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
        for(int idx=0;idx < nfp0;idx++)rho[idx] = ct.mix*rho[idx] + (1.0-ct.mix)*trho[idx];
    }
    else
    {
        if(ct.scf_steps <ct.freeze_orbital_step)
        {
            steps = ct.scf_steps - ct.freeze_rho_steps;
        }
        else
        {
            if(ct.charge_pulay_order ==1 )  ct.charge_pulay_order++;
            steps = ct.scf_steps - ct.freeze_orbital_step;
        }

        for (idx = 0; idx < nfp0; idx++)
        {
            tem = rho_old[idx];
            rho_old[idx] = rho[idx] - rho_old[idx];
            rho[idx] = tem;
        }

        if(ct.spin_flag)
            get_rho_oppo(rho, rho_oppo);
        //get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
        if(ct.scf_steps >= ct.freeze_rho_steps)
            Pulay_rho->Mixing(rho, rho_old);

    }
    delete(RT3);

    if(ct.spin_flag) get_rho_oppo(rho, rho_oppo);

    /* Update potential */
    RmgTimer *RT4 = new RmgTimer("2-SCF: update_pot");
    if(ct.scf_steps >= ct.freeze_rho_steps)
        UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);
    delete(RT4);

    CheckConvergence(vxc, vh, vxc_old, vh_old, rho, rho_pre, CONVERGENCE);

    RmgTimer *RT5 = new RmgTimer("2-SCF: get_te");
    delete(RT5);

    /* Update the orbitals */
    if(ct.scf_steps < ct.freeze_orbital_step)
    {
        steps = ct.scf_steps;
        if(ct.scf_steps == ct.freeze_rho_steps ) 
            ct.restart_mix = 1;
        else
            ct.restart_mix = 0;
        
        if(ct.LocalizedOrbitalLayout == LO_projection)
        {
            CalculateResidual(*LocalOrbital, *H_LocalOrbital, *LocalProj,  vtot_c, theta, Kbpsi_mat);
            H_LocalOrbital->WriteOrbitals(std::string(ct.outfile), Rmg_G);
            read_orbitals_on(ct.outfile, states1);
        }
        RmgTimer *RT6 = new RmgTimer("2-SCF: OrbitalOptimize");
        OrbitalOptimize(states, states1, vxc, vh, vnuc, rho, rhoc, vxc_old, vh_old);

        delete(RT6);
    }

    delete(RT);

    delete [] trho;
    delete [] rho_pre;
}                               /* end scf */



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

/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"
#include "common_prototypes1.h"
#include "Scalapack.h"
#include "prototypes_on.h"

//#include "main.h"
//#include "init_var.h"
//#include "RmgTimer.h"


#if !ELEMENTAL_LIBS
#include "blas.h"



void DiagScalapack(STATE *states, int numst, double *Hij_00, double *Bij_00, double *rho_matrix, double *theta_ptr)
{

    RmgTimer  *RT0 = new RmgTimer("3-DiagScalapack");
    int ione = 1;    /* blas constants */
    char *uplo = "u", *jobz = "v";

    int info;
    double zero = 0., one = 1., alpha;
    int st1, st_g;
    double *eigs;

    /* for parallel libraries */
    int mb= pct.desca[4];
    int  mxllda2;
    mxllda2 = MXLLDA * MXLCOL;
    double *work1 = new double[mxllda2]();



    RmgTimer *RT = new RmgTimer("3-DiagScalapack: pdsygvx ");
    /* If I'm in the process grid, execute the program */
    if (pct.scalapack_myrow < 0)
    {  
        printf("\n nprow, npcol %d %d", pct.scalapack_nprow, pct.scalapack_npcol);
        printf("\n we should use all proc for diag. somthing wrong");
        printf("\n gridpe = %d", pct.gridpe);
        exit(0);
    }


    /* 
     * SOLVE THE GENERALIZED EIGENVALUE PROBLEM:  m * z = lambda * matS * z 
     */

    /* Transform the generalized eigenvalue problem to a sStandard form */
    dcopy(&mxllda2, Hij, &ione, uu_dis, &ione);
    dcopy(&mxllda2, matB, &ione, l_s, &ione);
//dcopy(&mxllda2, matB, &ione, l_s, &ione);
//pdtran(&numst, &numst, &half, l_s, &ione, &ione, pct.desca,
//           &half, matB, &ione, &ione, pct.desca);

    char *range = "a";
    double vx = 0.0;
    double tol = 0.0;
    int eigs_found, eigvs_found;
    double orfac = 0.0;
    int *iwork, *ifail, *iclustr, lwork;
    double *gap, lwork_tmp, *work2;
    int liwork_tmp, liwork;

    ifail = new int[numst];
    eigs = new double[numst];
    iclustr = new int[2 * pct.scalapack_nprow * pct.scalapack_npcol];
    gap = new double[pct.scalapack_nprow * pct.scalapack_npcol];

    lwork = -1;
    liwork = -1;

    pdsygvx (&ione, jobz, range, uplo, &numst, uu_dis, &ione, &ione, pct.desca,
            l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
            &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
            &liwork_tmp, &liwork, ifail, iclustr, gap, &info);

    if (info)
    {
        printf ("\n pdsygvx query failed, info is %d", info);
        exit(0);
    }


    /*set lwork and liwork */
    lwork = (int) lwork_tmp + 1;
    liwork = liwork_tmp;

    work2 = new double[lwork];
    iwork = new int[liwork];




    tol = 1.0e-15;
    pdsygvx (&ione, jobz, range, uplo, &numst, uu_dis, &ione, &ione, pct.desca,
            l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
            &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, work2, &lwork,
            iwork, &liwork, ifail, iclustr, gap, &info);


    if (info)
    {
        printf ("\n pdsygvx failed, info is %d", info);
        exit(0);
    }



    delete(RT);


    RmgTimer *RT1 = new RmgTimer("3-DiagScalapack: calc_occ");
    for (st1 = 0; st1 < numst; st1++)
    {
        states[st1].eig[0] = eigs[st1];
    }

    delete [] ifail;
    delete [] iclustr;
    delete [] gap;
    delete [] work2;
    delete [] iwork;
    delete [] eigs;

    delete(RT1);

    RmgTimer *RT5 = new RmgTimer("3-DiagScalapack: pscal occ ");
    //   uu_dis = zz_dis *(occ_diag)
    dcopy(&mxllda2, zz_dis, &ione, uu_dis, &ione);
    for(st1 = 0; st1 <  MXLCOL; st1++)
    {

        st_g = (st1/mb) * pct.scalapack_npcol * mb + pct.scalapack_mycol *mb +st1%mb;

        if(st_g >= numst) 
            alpha = 0.0;
        else
            alpha = states[st_g].occupation[0];
        dscal(&MXLLDA, &alpha, &uu_dis[st1 * MXLLDA], &ione);
    }

    delete(RT5);

    RmgTimer *RT3 = new RmgTimer("3-DiagScalapack: gemm ");

    pdgemm("N", "T", &numst, &numst, &numst, &one,
            uu_dis, &ione, &ione, pct.desca,
            zz_dis, &ione, &ione, pct.desca, &zero, mat_X, &ione, &ione, pct.desca);

    delete(RT3);


    RmgTimer *RT1b = new RmgTimer("3-DiagScalapack: (S^-1)H");

    int *ipiv = new int[numst];
    /* Compute matrix theta = matB^-1 * Hij  */
    pdgetrf(&numst, &numst, matB, &ione, &ione, pct.desca, ipiv, &info);
    if(info !=0)
    { 
        printf("\n error in pdgetrf in mg_eig.c INFO = %d\n", info);
        fflush(NULL);
        exit(0);
    }

    dcopy(&mxllda2, Hij, &ione, uu_dis, &ione);
    pdgetrs("N", &numst, &numst, matB, &ione, &ione, pct.desca, ipiv, 
            uu_dis, &ione, &ione, pct.desca, &info);

    delete [] ipiv;


    double t1 = 2.0;
    dscal(&mxllda2, &t1, uu_dis, &ione);

  //  Cpdgemr2d(numst, numst, uu_dis, ione, ione, pct.desca, theta_ptr, ione, ione,
  //          pct.descb, pct.desca[1]);
    delete(RT1b);

    delete [] work1;


    delete(RT0);
}
#endif

