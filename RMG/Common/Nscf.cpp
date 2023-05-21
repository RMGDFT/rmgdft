/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "Solvers.h"
#include "RmgParallelFft.h"
#include "RmgSumAll.h"




template bool Nscf<double> (double *, double *, double *, double *,
          double *, double *, double *, double *, double *, double *, int ,
          int , Kpoint<double> **, std::vector<double> &);
template bool Nscf<std::complex<double> > (double *, double *, double *, double *,
          double *, double *, double *, double *, double *, double *, int ,
          int , Kpoint<std::complex<double>> **, std::vector<double> &);

template <typename OrbitalType> bool Nscf (double * vxc, double *vxc_in, double * vh, double *vh_in, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV)
{

    RmgTimer RT0("2-Scf steps");
    RmgTimer RTt("1-TOTAL: run: Scf steps");

    bool CONVERGED = false;
    static double eigsum_old = 0.0;
    double *vtot, *vtot_psi;
    double *vxc_psi=NULL;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;
    int P0_BASIS;

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    /* allocate memory for eigenvalue send array and receive array */
    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];

    /* save old vhxc + vnuc */
    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
    }

    if(ct.noncoll) vxc_psi = new double[4*P0_BASIS];

    // Transfer vtot from the fine grid to the wavefunction grid
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
    if(ct.noncoll)
        for(int is = 0; is < ct.nspin; is++)
        {
            GetVtotPsi (&vxc_psi[is*P0_BASIS], &vxc[is*FP0_BASIS], Rmg_G->default_FG_RATIO);
        }

    /*Generate the Dnm_I */
    get_ddd (vtot, vxc, true);


    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {

        //if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) || 
        //        ((ct.scf_steps < ct.davidson_premg) && (ct.md_steps == 0) && (ct.runflag != RESTART )) ||
        //        (ct.xc_is_hybrid && Functional::is_exx_active())) 
        {
            RmgTimer *RT1 = new RmgTimer("2-Scf steps: MgridSubspace");
            Kptr[kpt]->MgridSubspace(vtot_psi, vxc_psi);
            delete RT1;
        }
       // else if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) {
       //     int notconv;
       //     RmgTimer *RT1 = new RmgTimer("2-Scf steps: Davidson");
       //     Kptr[kpt]->Davidson(vtot_psi, vxc_psi, notconv);
       //     delete RT1;
       // }

        // Needed to ensure consistency with some types of kpoint parrelization
        MPI_Barrier(pct.grid_comm);
    } // end loop over kpoints

    if (ct.nspin == 2)
        GetOppositeEigvals (Kptr);

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {
        int factor = 1;
        if(ct.nspin == 2) factor = 2;
        ct.kp[kpt+pct.kstart].eigs.resize(ct.num_states * factor);
        for(int st = 0; st < ct.num_states; st++)
            ct.kp[kpt+pct.kstart].eigs[st] = Kptr[kpt]->Kstates[st].eig[0];

        if(ct.nspin == 2)
        {
            for(int st = 0; st < ct.num_states; st++)
                ct.kp[kpt+pct.kstart].eigs[st + ct.num_states] = Kptr[kpt]->Kstates[st].eig[1];
        }
    }


    if(ct.ldaU_mode != LDA_PLUS_U_NONE && (ct.ns_occ_rms >1.0e-15 || ct.scf_steps ==0) )
    {
        RmgTimer("3-MgridSubspace: ldaUop x psi");
        int pstride = Kptr[0]->ldaU->ldaU_m;
        int num_ldaU_ions = Kptr[0]->ldaU->num_ldaU_ions;
        int occ_size = ct.nspin * num_ldaU_ions * pstride * pstride;

        doubleC_4d_array old_ns_occ, new_ns_occ;
        old_ns_occ.resize(boost::extents[ct.nspin][num_ldaU_ions][Kptr[0]->ldaU->ldaU_m][Kptr[0]->ldaU->ldaU_m]);
        new_ns_occ.resize(boost::extents[ct.nspin][num_ldaU_ions][Kptr[0]->ldaU->ldaU_m][Kptr[0]->ldaU->ldaU_m]);

        old_ns_occ =  Kptr[0]->ldaU->ns_occ;

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            LdaplusUxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->orbitalsint_local);
            Kptr[kpt]->ldaU->calc_ns_occ(Kptr[kpt]->orbitalsint_local, 0, Kptr[kpt]->nstates);
            for(int idx = 0; idx< occ_size; idx++)
            {
                new_ns_occ.data()[idx] += Kptr[kpt]->ldaU->ns_occ.data()[idx];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)new_ns_occ.data(), occ_size * 2, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

        Rmg_Symm->symm_nsocc(new_ns_occ.data(), pstride, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
        if(ct.AFM)
        {
            Rmg_Symm->nsocc_AFM(new_ns_occ, Kptr[0]->ldaU->ldaU_m, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
        }


        // for spin-polarized case, only the first spin on the pe are symmetrized so communicate to opposite spin
        if(ct.nspin == 2 && !ct.AFM)
        {
            MPI_Status status;
            int len = 2 * Kptr[0]->ldaU->num_ldaU_ions * pstride * pstride;
            double *sendbuf = (double *)new_ns_occ.data();
            double *recvbuf = sendbuf + len;
            MPI_Sendrecv(sendbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
                    recvbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);

        }


        ct.ns_occ_rms = 0.0;
        for(int idx = 0; idx < occ_size; idx++)
        {
            ct.ns_occ_rms += std::norm (old_ns_occ.data()[idx] - new_ns_occ.data()[idx]);
        }
        ct.ns_occ_rms = sqrt(ct.ns_occ_rms / occ_size);

        MixLdaU(occ_size *2, (double *)new_ns_occ.data(), (double *)old_ns_occ.data(), Kptr[0]->ControlMap, false);

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            Kptr[kpt]->ldaU->ns_occ = old_ns_occ;
            Kptr[kpt]->ldaU->calc_energy();
        }
        if(ct.verbose) Kptr[0]->ldaU->write_ldaU();

    }

    double eigsum = 0.0;

    int nspin = 1;
    if(ct.nspin == 2 && !ct.AFM) nspin = 2;
    for (int idx = 0; idx < nspin; idx++)
    {
        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {

            double t1 = 0.0;
            for (int state = 0; state < ct.num_states; state++)
            {

                t1 += (Kptr[kpt]->Kstates[state].occupation[idx] *
                        Kptr[kpt]->Kstates[state].eig[idx]);

            }
            eigsum += t1 * Kptr[kpt]->kp.kweight;
        }
    }

    eigsum = RmgSumAll(eigsum, pct.kpsub_comm);

    ct.scf_accuracy = eigsum - eigsum_old;
    eigsum_old = eigsum;

    ct.rms = ct.scf_accuracy;
    if(fabs(ct.scf_accuracy) <  ct.adaptive_thr_energy) CONVERGED = true;


    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
#if !(__CYGWIN__ || __MINGW64__ || _WIN32 || _WIN64)
    fsync( fileno(ct.logfile) );
#endif

    /* release memory */
    delete [] vtot;
    delete [] vtot_psi;
    if(ct.noncoll) delete [] vxc_psi;

    return CONVERGED;
}                               /* end scf */


/******/
