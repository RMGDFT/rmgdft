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
#include "LdaU_on.h"
#include "GpuAlloc.h"
#include "Exx_on.h"
#include "BlockTriMatrix.h"
#include "blas_driver.h"
#include "Symmetry.h"

#define DELTA_V_MAX 1.0

void update_pot(double *, double *, double *, double *, double *, double *, double *,
        double *, double *, int *, STATE * states);


void Scf_on_proj(STATE * states, double *vxc, double *vh,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore,
        double * vxc_old, double * vh_old, int *CONVERGENCE, bool freeze_orbital,
        double *rho_pre, double *trho, double *rho_matrix_local, double *theta_local,
        double *CC_res_local, double *Hij_local, double *Sij_local,
        double *Hij_glob, double *Sij_glob)
{
    int idx, ione = 1;
    double tem;
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int pbasis= Rmg_G->get_P0_BASIS(1);

    RmgTimer *RT0;


    RmgTimer *RT = new RmgTimer("2-SCF");

    MemcpyHostDevice(LocalOrbital->storage_size, LocalOrbital->storage_cpu, LocalOrbital->storage_gpu);

    RT0 = new RmgTimer("2-SCF: Kbpsi");


    RmgTimer *RT1 = new RmgTimer("2-SCF: Kbpsi: LOxLO");
    LO_x_LO(*LocalProj, *LocalOrbital, Kbpsi_mat_local, *Rmg_G);
    delete RT1;


    MPI_Barrier(MPI_COMM_WORLD);
    //  now Kbpsi_mat_local stores the results from different processors, need to be summer over all processors.
    RT1 = new RmgTimer("2-SCF: Kbpsi: local_to_global");
    mat_local_to_glob(Kbpsi_mat_local, Kbpsi_mat, *LocalProj, *LocalOrbital, 0, LocalProj->num_tot, 
            0, LocalOrbital->num_tot, true);
    delete RT1;

    MPI_Barrier(MPI_COMM_WORLD);
    RT1 = new RmgTimer("2-SCF: Kbpsi: local_to_tri");
    //BlockTriMatrix<double> Kbpsi_tri(ct.num_blocks, LocalProj->num_tot, LocalOrbital->num_tot, ct.block_dim_nl, ct.block_dim_phi, true);

    //Kbpsi_tri.Local2BlockTri(Kbpsi_mat_local,  *LocalProj, *LocalOrbital);
    delete RT1;

    RT1 = new RmgTimer("2-SCF: Kbpsi: global_to_local");
    mat_global_to_local(*LocalProj, *LocalOrbital, Kbpsi_mat, Kbpsi_mat_local); 
    delete RT1;

    // now Kbpsi_mat_local store the correct values.
    delete RT0;

    RT0 = new RmgTimer("2-SCF: HS mat");
    // mat_X charge density matrix in distributed way
    // uu_dis theta = (S^-1 H) in distributed way.

    int num_orb = LocalOrbital->num_thispe;

    int num_tot = LocalOrbital->num_tot;
    RT1 = new RmgTimer("2-SCF: HS mat: ApplyH");
    ApplyHphi(*LocalOrbital, *H_LocalOrbital, vtot_c);
    MemcpyHostDevice(H_LocalOrbital->storage_size, H_LocalOrbital->storage_cpu, H_LocalOrbital->storage_gpu);
    delete RT1;

    RT1 = new RmgTimer("2-SCF: HS mat: LOxLO");
    LO_x_LO(*LocalOrbital, *H_LocalOrbital, Hij_local, *Rmg_G);
    LO_x_LO(*LocalOrbital, *LocalOrbital, Sij_local, *Rmg_G);

    delete RT1;

    RT1 = new RmgTimer("2-SCF: HS mat: reduce");
    mat_local_to_glob(Hij_local, Hij_glob, *LocalOrbital, *LocalOrbital, 0, num_tot, 0, num_tot, false);
    mat_local_to_glob(Sij_local, Sij_glob, *LocalOrbital, *LocalOrbital, 0, num_tot, 0, num_tot, false);
    delete RT1;

    RT1 = new RmgTimer("2-SCF: HS mat: Hvnl");
    GetHvnlij_proj(Hij_glob, Sij_glob, Kbpsi_mat, Kbpsi_mat, 
            num_tot, num_tot, LocalProj->num_tot, true);
    delete RT1;

    RT1 = new RmgTimer("2-SCF: HS mat: reduce");
    int sum_dim = num_tot * num_tot;
    MPI_Allreduce(MPI_IN_PLACE, Hij_glob, sum_dim, MPI_DOUBLE, MPI_SUM, LocalOrbital->comm);
    MPI_Allreduce(MPI_IN_PLACE, Sij_glob, sum_dim, MPI_DOUBLE, MPI_SUM, LocalOrbital->comm);


    delete RT1;

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_glob, 6, num_tot);
        print_matrix(Sij_glob, 6, num_tot);

    }
    if(ct.xc_is_hybrid && Functional::is_exx_active())
    {
         Exx_onscf->HijExx(Hij_glob, *LocalOrbital);
    }
    delete RT0;


    switch(ct.subdiag_driver) {
        case SUBDIAG_CUSOLVER:
            {
                RT0 = new RmgTimer("2-SCF: DiagGpu");
                DiagGpu(states, ct.num_states, Hij_glob, Sij_glob, rho_matrix_local, theta_local, CC_res_local);
                delete RT0;
                break;
            }
        default:
            {

                RT0 = new RmgTimer("2-SCF: DiagScalapack");
                for(int i = 0; i < num_orb; i++) 
                {
                    for(int j = 0; j < num_orb; j++) 
                    {
                        if(i == j ) CC_res_local[i*num_orb + j] = 1.0;
                        else CC_res_local[i*num_orb + j] = 0.0;
                    }
                }
                mat_global_to_dist(Hij, pct.desca, Hij_glob);
                mat_global_to_dist(matB, pct.desca, Sij_glob);

                if(ct.is_gamma)
                    DiagScalapack<double>(states, ct.num_states, Hij, matB);
                else
                    DiagScalapack<std::complex<double>>(states, ct.num_states, Hij, matB);


                mat_dist_to_local(mat_X, pct.desca, rho_matrix_local, *LocalOrbital);
                mat_dist_to_local(uu_dis, pct.desca, theta_local, *LocalOrbital);
                delete RT0;
            }
    }


    if(ct.num_ldaU_ions > 0)
    {
        RT0 = new RmgTimer("2-SCF: LDA+U");
        ldaU_on->calc_ns_occ(*LocalOrbital, mat_X, *Rmg_G);
        if(ct.verbose) ldaU_on->write_ldaU();
        delete RT0;

    }

    get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);

    if (pct.gridpe == 0 && ct.occ_flag == 1)
        rmg_printf("FERMI ENERGY = %15.8f\n", ct.efermi * Ha_eV);

    dcopy(&nfp0, rho, &ione, rho_old, &ione);
    dcopy(&nfp0, rho, &ione, rho_pre, &ione);

    RT0 = new RmgTimer("2-SCF: GetNewRho");
    for (idx = 0; idx < nfp0; idx++)trho[idx] = rho[idx];

    if(ct.scf_steps >= ct.freeze_rho_steps)
    {
        GetNewRho_proj(*LocalOrbital, *H_LocalOrbital, rho, rho_matrix_local);
    }
    if(ct.spin_flag)
        get_rho_oppo(rho, rho_oppo);
    if (ct.AFM)
    {
        RmgTimer RTT("symm_rho");
        Rmg_Symm->symmetrize_rho_AFM(rho, rho_oppo);
    }
    delete RT0;

    RT0 = new RmgTimer("2-SCF: Rho mixing");
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
        printf("\n Warning: total charge Normalization constant = %15.12e  \n", t2-1.0);


    if(ct.charge_mixing_type == 0)
    {
        if(ct.spin_flag)
            get_rho_oppo(rho, rho_oppo);
        //     get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
        for(int idx=0;idx < nfp0;idx++)rho[idx] = ct.mix*rho[idx] + (1.0-ct.mix)*trho[idx];
    }
    else
    {

        for (idx = 0; idx < nfp0; idx++)
        {
            tem = rho_old[idx];
            rho_old[idx] = rho[idx] - rho_old[idx];
            rho[idx] = tem;
        }

        if(ct.spin_flag)
            get_rho_oppo(rho, rho_oppo);
        //get_te(rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);
        if(ct.scf_steps ==ct.freeze_orbital_step)
        {
            if(ct.charge_pulay_order ==1 )  ct.charge_pulay_order++;
            Pulay_rho->Refresh();
        }

        if(ct.scf_steps >= ct.freeze_rho_steps)
        {
            Pulay_rho->Mixing(rho, rho_old);
        }

    }

    delete RT0;

    if(ct.spin_flag) get_rho_oppo(rho, rho_oppo);

    /* Update potential */
    RmgTimer *RT4 = new RmgTimer("2-SCF: update_pot");
    if(ct.scf_steps >= ct.freeze_rho_steps)
        UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);
    delete(RT4);

    CheckConvergence(vxc, vh, vxc_old, vh_old, rho, rho_pre, CONVERGENCE);

    /* Update the orbitals */
    if(!freeze_orbital && ct.is_gamma)
    {
        if(ct.scf_steps == ct.freeze_rho_steps ) 
            ct.restart_mix = 1;
        else
            ct.restart_mix = 0;

        RT0 = new RmgTimer("2-SCF: Residual calculation");
        CalculateResidual(*LocalOrbital, *H_LocalOrbital, *LocalProj,  vtot_c, theta_local, Kbpsi_mat, CC_res_local);
        delete RT0;
        my_sync_device();
        for(int st = 0; st < LocalOrbital->num_thispe; st++)
        {

            LocalOrbital->ApplyBoundary(&H_LocalOrbital->storage_cpu[st * pbasis], st);
            //for(int idx = 0; idx < pbasis; idx++) if (!LocalOrbital->mask[st * pbasis + idx])
            //    H_LocalOrbital->storage_cpu[st * pbasis + idx] = 0.0;
        }
        RT0 = new RmgTimer("2-SCF: orbital precond and mixing");

        Pulay_orbital->Mixing(LocalOrbital->storage_cpu, H_LocalOrbital->storage_cpu);
        //Pulay_orbital->MixingOrbitalBroyden(LocalOrbital->storage_cpu, H_LocalOrbital->storage_cpu, vh_out, vh_in);
        RmgTimer *RT1 = new RmgTimer("2-SCF: orbital precond and mixing: normalize");
        LocalOrbital->Normalize();
        delete RT1;
        delete RT0;

    }
    delete(RT);

}                               /* end scf */



