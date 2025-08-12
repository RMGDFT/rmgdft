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

#include <complex>
#include <boost/pool/pool.hpp>
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "RmgSumAll.h"
#include "BlasWrappers.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "packfuncs.h"
#include "transition.h"
#include "BaseThread.h"
#include "MpiQueue.h"
#include "GatherScatter.h"
#include "Solvers.h"
#include "rmg_complex.h"



template <typename T>
double ComputeEig(int n, T *A, T *B, T *D)
{
    double s1[2];
    s1[0] = 0.0;
    s1[1] = 0.0;

    for(int idx = 0;idx < n;idx++)
    {
        s1[0] += std::real(A[idx])*std::real(B[idx]) + std::imag(A[idx])*std::imag(B[idx]);
        s1[1] += std::real(A[idx])*std::real(D[idx]) + std::imag(A[idx])*std::imag(D[idx]);
    }

    int length = 2;
    GlobalSums (s1, length, pct.coalesced_grid_comm);
    return  s1[0] / s1[1];

}

std::atomic<bool> reduce_it = false;

template void MgEigState<double,float>(Kpoint<double> *, State<double> *, double *, double *, double *, double *, double *, int);
template void MgEigState<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, double *, double *, int);
template void MgEigState<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *, double *,
        double *, std::complex<double> *, std::complex<double> *, int);
template void MgEigState<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *, double *,
        double *, std::complex<double> *, std::complex<double> *, int);


template <typename OrbitalType, typename CalcType>
void MgEigState (Kpoint<OrbitalType> *kptr, State<OrbitalType> * sp, double * vtot_psi, double *coarse_vtot, double *vxc_psi, OrbitalType *nv, OrbitalType *ns, int vcycle)
{
    RmgTimer RT("Mg_eig");
    BaseThread *Thread = BaseThread::getBaseThread(0);
    int tid = Thread->get_thread_tid();

    bool freeze_occupied = true;

    // We want a clean exit if user terminates early
    CheckShutdown();

    // We can't just skip the occupied orbitals if they are frozen since we process states in blocks and
    // combine communications. So for now set a flag indicating whether we update the orbital or not. It should
    // be possible to fix this at a higher level at some point though so unneccessary work is not done.
    if(Verify ("freeze_occupied", true, kptr->ControlMap) && (sp->occupation[0] > 0.0)) freeze_occupied = false;

    if(Verify ("calculation_mode", "Band Structure Only", kptr->ControlMap) )
        freeze_occupied = true;

    bool using_davidson = Verify ("kohn_sham_solver","davidson", kptr->ControlMap);

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;
    TradeImages *T = kptr->T;

    double eig=0.0, diag, t1;
    int eig_pre[MAX_MG_LEVELS] = { 0, 8, 8, 8, 8, 8, 8, 8 };
    int eig_post[MAX_MG_LEVELS] = { 0, 8, 8, 8, 8, 8, 8, 8 };
    int potential_acceleration;
    Mgrid MG(L, T);

    // Once the SCF error gets small increase this to remove
    // interpolation errors from coarse to fine.
    int gl_pst = ct.eig_parm.gl_pst;
    if(ct.scf_accuracy < 1.0e-9) gl_pst=3;

    int nits = ct.eig_parm.gl_pre + gl_pst;
    int dimx = G->get_PX0_GRID(1) * pct.coalesce_factor;
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int NX_GRID = G->get_NX_GRID(1);
    int NY_GRID = G->get_NY_GRID(1);
    int NZ_GRID = G->get_NZ_GRID(1);

    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);
    int levels = ct.eig_parm.levels;
    bool do_mgrid = true;
    double mg_step = ct.eig_parm.sb_step;
    double fg_step = ct.eig_parm.mg_timestep;

    if(reduce_it) fg_step = std::min(2.0/3.0, ct.eig_parm.mg_timestep);

    double gl_step = ct.eig_parm.gl_step;
    if ((ct.runflag == RANDOM_START) && (ct.scf_steps < 2)) do_mgrid = false;

    double Zfac = 2.0 * ct.max_zvalence;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);
    int pbasis_noncoll = pbasis * ct.noncoll_factor;

    size_t aratio = sizeof(OrbitalType) / sizeof(CalcType);

    // Per thread memory pool in Kpoint class. If you need more than 32 chunks modify
    // pool creation code in Kpoint.cpp. Done like this to enable changing the type
    // of memory allocated in one place.
    int pool_blocks = 0;
    boost::pool<rmg_user_allocator> *p = kptr->kalloc[tid];
    CalcType *res2_t = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *work2_t = (CalcType *)p->ordered_malloc(4);pool_blocks+=4;
    CalcType *work1_t = (CalcType *)p->ordered_malloc(4);pool_blocks+=4;
    CalcType *sg_twovpsi_t  =  (CalcType *)p->ordered_malloc(4);pool_blocks+=4;
    OrbitalType *saved_psi  = (OrbitalType *)p->ordered_malloc(aratio);pool_blocks+=aratio;
    double *dvtot_psi = (double *)p->ordered_malloc(aratio);pool_blocks+=aratio;
    double *c_vtot = (double *)p->ordered_malloc(2);pool_blocks+=2;
    CalcType *tmp_psi_t  = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *res_t  =  (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *twork_t  = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    OrbitalType *nv_t  = (OrbitalType *)p->ordered_malloc(aratio);pool_blocks+=aratio;

    // Copy double precision psi into correct precison array
    GatherPsi(G, pbasis_noncoll, sp->istate, kptr->orbital_storage, tmp_psi_t, pct.coalesce_factor);

    // Copy nv into local array
    if(ct.coalesce_states)
        GatherPsi(G, pbasis_noncoll, sp->istate, nv, nv_t, pct.coalesce_factor);
//        GatherPsi(G, pbasis_noncoll, sp->istate, kptr->nv, nv_t, pct.coalesce_factor);
    else
        GatherPsi(G, pbasis_noncoll, 0, nv, nv_t, 1);

    // Set up coarse vtot
    for(int idx=0;idx < pbasis;idx++) c_vtot[idx] = -coarse_vtot[idx];

    // For USPP copy double precision ns into correct precision temp array. For NCPP ns=psi. */
    if(ct.norm_conserving_pp)
    {
        for(int ix=0;ix < pbasis_noncoll;ix++) work1_t[ix] = tmp_psi_t[ix];
    }
    else
    {
        GatherPsi(G, pbasis_noncoll, sp->istate, kptr->ns, work1_t, pct.coalesce_factor);
    }

    /* Save in res2 */
    for(int ix=0;ix < pbasis_noncoll;ix++) res2_t[ix] = work1_t[ix];

    // Setup some potential acceleration stuff
        for(int idx = 0;idx <pbasis_noncoll;idx++) saved_psi[idx] = (OrbitalType)tmp_psi_t[idx];
    potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    if(potential_acceleration) {
        PotentialAccelerationWait(sp->istate, kptr->nstates, ct.dvh_skip);
    }

    /* Smoothing cycles */
    for (int cycles = 0; cycles <= nits; cycles++)
    {
        /* Apply Hamiltonian */
        {
            RmgTimer RT1("Mg_eig: apply hamiltonian");
            diag=ApplyHamiltonian<OrbitalType,CalcType> (kptr, sp, sp->istate, tmp_psi_t, work1_t, vtot_psi, vxc_psi, nv_t, false);
        }

        // Copy saved application to ns to res
        memcpy(res_t, res2_t, pbasis_noncoll * sizeof(CalcType));

        /* If this is the first time through compute the eigenvalue */
        if ((cycles == 0) || (potential_acceleration != 0) || (using_davidson && (cycles == 0))) 
        {
            eig = ComputeEig(pbasis_noncoll, tmp_psi_t, work1_t, res_t);
            // Save this for variational energy correction
            if((cycles == 0) && (vcycle == 0)) sp->feig[0]=eig;


            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
             * are not defined yet*/
            if(freeze_occupied) {

                if ((ct.diag == 1) && (potential_acceleration == 0))
                {
                    if ((ct.scf_steps == 0) && (ct.exx_steps == 0))
                    {
                        sp->eig[0] = eig;
                        sp->oldeig[0] = eig;
                    }
                    else
                        eig = sp->eig[0];
                }
                else
                {
                    sp->eig[0] = eig;
                    if(ct.scf_steps == 0 && (ct.exx_steps == 0)) {
                        sp->oldeig[0] = eig;
                    }
                }

                if(potential_acceleration) {
                    t1 = eig;
                    double s1 = std::pow(0.5, cycles+1);
                    eig = s1 * eig + (1-s1) * sp->oldeig[0];
                    sp->oldeig[0] = t1;
                }
            }
            else {
                eig = sp->eig[0];
            }

        }
        // Get the residual
        CalcType f1(2.0*eig);
        for (int idx = 0; idx <pbasis_noncoll; idx++) res_t[idx] = f1 * res_t[idx] - 2.0*work1_t[idx];
        sp->res[cycles] = 0.0;
        for (int idx = 0; idx <pbasis_noncoll; idx++) sp->res[cycles] += std::real(res_t[idx] * std::conj(res_t[idx]));
        GlobalSums (&sp->res[cycles], 1, pct.coalesced_grid_comm);
        sp->res[cycles] *= get_vel();

        if(cycles > 0 && sp->res[cycles] > sp->res[cycles-1] && cycles < ct.eig_parm.gl_pre)
        {
            if(ct.verbose && pct.gridpe==0)
                printf("\nResidual increasing  %d   %12.8e  %12.8e.\n",sp->istate, sp->res[cycles-1], sp->res[cycles]);
            fg_step = std::min(2.0/3.0, ct.eig_parm.mg_timestep);
            reduce_it = true;
            if(ct.verbose && pct.gridpe==0) printf("\nNot smoothing! Resetting state %d\n", sp->istate);
            for(int idx = 0;idx <pbasis_noncoll;idx++) tmp_psi_t[idx] = saved_psi[idx];
            break;
        }
        for(int is = 0; is < ct.noncoll_factor; is++)
        {
            /* Now either smooth the wavefunction or do a multigrid cycle */
            if ((cycles == ct.eig_parm.gl_pre) && do_mgrid )
            {

                /* Do multigrid step with solution returned in sg_twovpsi */
                {

                    RmgTimer RT1("Mg_eig: mgrid_solv");

                    // We use a residual correction multigrid scheme where the right hand side is the residual
                    // so single precision is adequate for the correction since the errors from lower precision
                    // will be approximately 7 decimal digits smaller than the original error we are correcting
                    // for. The std::conditional_t typdefs allow us to cleanly handle the multiple combinations
                    // of OrbitalType and CalcType. The combinations are as follows.
                    //
                    //   OrbitalType = double, CalcType = double
                    //   OrbitalType = double, CalcType = float
                    //   OrbitalType = std::complex<double>, CalcType = std::complex<double>
                    //   OrbitalType = std::complex<double>, CalcType = std::complex<float>
                    // with mg_type always float or std::complex<float>
                    typedef typename std::conditional_t< std::is_same<CalcType, double>::value, float,
                                     std::conditional_t< std::is_same<CalcType, std::complex<double>>::value, std::complex<float>,
                                     std::conditional_t< std::is_same<CalcType, std::complex<float>>::value, std::complex<float>, float> > > mgtype_t;
                    typedef typename std::conditional_t< std::is_same<CalcType, double>::value, double,
                                     std::conditional_t< std::is_same<CalcType, std::complex<double>>::value, std::complex<double>,
                                     std::conditional_t< std::is_same<CalcType, std::complex<float>>::value, std::complex<float>, float> > > convert_type_t;

                    mgtype_t *v_mat = (mgtype_t *)&sg_twovpsi_t[sbasis];
                    mgtype_t *f_mat = (mgtype_t *)&work1_t[sbasis];
                    mgtype_t *twork_tf = (mgtype_t *)twork_t;
                    mgtype_t *work2_tf = (mgtype_t *)work2_t;


                    int ixoff, iyoff, izoff;
                    int dx2 = MG.MG_SIZE (dimx, 0, NX_GRID, G->get_PX_OFFSET(1), dimx, &ixoff, ct.boundaryflag);
                    int dy2 = MG.MG_SIZE (dimy, 0, NY_GRID, G->get_PY_OFFSET(1), dimy, &iyoff, ct.boundaryflag);
                    int dz2 = MG.MG_SIZE (dimz, 0, NZ_GRID, G->get_PZ_OFFSET(1), dimz, &izoff, ct.boundaryflag);

                    if((dx2 < 0) || (dy2 < 0) || (dz2 < 0)) {
                        printf("Multigrid error: Grid cannot be coarsened. Most likely the current grid is not divisable by 2 or 4. It is recommended to use grid that is, at minimum, divisable by 4. The current grid is %d %d %d" , NX_GRID, NY_GRID, NZ_GRID);
                        exit(0);
                    }

                    /* Pack the residual data into multigrid array */
                    CPP_pack_ptos_convert (twork_tf, &res_t[is*pbasis], dimx, dimy, dimz);
                    T->trade_images (twork_tf, dimx, dimy, dimz, FULL_TRADE);
                    MG.mg_restrict (twork_tf, f_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

                    MG.mgrid_solv (v_mat, f_mat, work2_tf,
                            dx2, dy2, dz2, 2.0*hxgrid, 2.0*hygrid, 2.0*hzgrid, 
                            1, levels, eig_pre, eig_post, 1, 
                            mg_step, 2.0*Zfac, 0.0, NULL,
                            //ct.eig_parm.sb_step, Zfac, 0.0, c_vtot,
                            NX_GRID, NY_GRID, NZ_GRID,
                            G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                            dimx, dimy, dimz, ct.boundaryflag);

                    MG.mg_prolong (twork_tf, v_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
                    CopyAndConvert(sbasis, (mgtype_t *)twork_tf, (convert_type_t *)sg_twovpsi_t);

                }

                /* The correction is in a smoothing grid so we use this
                 * routine to update the orbital which is stored in a physical grid.
                 */

                t1 = -fg_step;
                CPP_pack_stop_axpy<CalcType> (sg_twovpsi_t, &tmp_psi_t[is*pbasis], t1, dimx, dimy, dimz);

            }
            else
            {
                for (int idx = 0; idx <pbasis; idx++)
                {
                    OrbitalType t5 = gl_step / std::abs(diag + vtot_psi[idx] + 0.5*kptr->kp.kmag);
                    tmp_psi_t[idx + is * pbasis] += t5 * (OrbitalType)res_t[idx + is * pbasis];
                }

                if (cycles == 0)
                {

                    // If occupied orbitals are frozen we compute residuals 
                    if(freeze_occupied) {
                        //double t2 = ZERO;
                        //for (int idx = 0; idx <pbasis; idx++) t2 += std::norm(res_t[idx]);
                        //GlobalSums (&t2, 1, pct.coalesced_grid_comm);
                        //t2 = RmgSumAll (t2, pct.coalesced_grid_comm);
                        //t1 = (double) (ct.psi_nbasis);
                        //sp->res = sqrt (t2 / t1);
                        //                    if(pct.imgpe == 0) std::cout << "Orbital " << sp->istate << " residual = " << sp->res << std::endl;
                    }

                }

            }
        }

    }                           /* end for */

    if(potential_acceleration)
        PotentialAcceleration(kptr, sp, vtot_psi, dvtot_psi, tmp_psi_t, saved_psi);

    // Copy single precision orbital back to double precision
    if(freeze_occupied)
        ScatterPsi(G, pbasis_noncoll, sp->istate, tmp_psi_t, kptr->orbital_storage, pct.coalesce_factor);

    p->free(res2_t, pool_blocks);

} // end MgEigState


