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

template double ComputeEig<double>(int n, double *A, double *B, double *D);
template double ComputeEig<std::complex<double>>(int n, std::complex<double> *A, std::complex<double> *B, std::complex<double> *D);
template double ComputeEig<float>(int n, float *A, float *B, float *D);
template double ComputeEig<std::complex<float>>(int n, std::complex<float> *A, std::complex<float> *B, std::complex<float> *D);

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
    BaseThread *Thread = BaseThread::getBaseThread(0);
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;
    int tid = Thread->get_thread_tid();

    // Save in case needed for variational energy correction term
    sp->feig[0]=sp->eig[0];

    // We want a clean exit if user terminates early
    CheckShutdown();

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;
    TradeImages *T = kptr->T;

    double eig=0.0, t1;
    int eig_pre[MAX_MG_LEVELS] = { 0, 6, 6, 6, 6, 6, 6, 6 };
    int eig_post[MAX_MG_LEVELS] = { 0, 6, 6, 6, 6, 6, 6, 6 };
    int potential_acceleration;
    Mgrid MG(L, T);

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

    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
    double diag = FD.fd_coeff0(ct.kohn_sham_fd_order, hxgrid);

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
    CalcType *tmp_psi_t  = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *res_t  =  (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *ihu_t  =  (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *r0_t  =  (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *hr0_t  =  (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *twork_t  = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    double *dinv  = (double *)p->ordered_malloc(1);pool_blocks++;
    CalcType *rmmres_t = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    OrbitalType *nv_t  = (OrbitalType *)p->ordered_malloc(aratio);pool_blocks+=aratio;

    // Copy double precision psi into correct precison array
    GatherPsi(G, pbasis_noncoll, sp->istate, kptr->orbital_storage, tmp_psi_t, pct.coalesce_factor);

    // Copy nv into local array
    if(ct.coalesce_states)
        GatherPsi(G, pbasis_noncoll, sp->istate, nv, nv_t, pct.coalesce_factor);
    else
        GatherPsi(G, pbasis_noncoll, 0, nv, nv_t, 1);

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
    for(int idx = 0;idx <pbasis_noncoll;idx++)
        saved_psi[idx] = (OrbitalType)tmp_psi_t[idx];

    potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    if(potential_acceleration) {
        PotentialAccelerationWait(sp->istate, kptr->nstates, ct.dvh_skip);
    }

    for(int is = 0; is < ct.noncoll_factor; is++)
    {
        for(int ix=0;ix < pbasis;ix++)
        {
            dinv[ix+is*pbasis] = 1.0 / std::abs(diag + 2.0*vtot_psi[ix] + 0.5*kptr->kp.kmag);
        }
    }

    bool is_jacobi = true;
    if(ct.norm_conserving_pp) is_jacobi = false;
    mgsmoother<OrbitalType,CalcType>(kptr, sp,
              tmp_psi_t, work1_t, res_t, ihu_t, r0_t,
              vtot_psi, vxc_psi, dinv,
              nv_t, res2_t,
              sp->eig[0], 3, is_jacobi, ct.lambda_max, ct.lambda_min, vcycle);

    if(ct.use_rmm_diis)
        sp->dptr->addfunc(saved_psi);

    // Check if residuals were decreasing and if not abort smoothing for
    // this state.
    bool smooth_status = (sp->res[0] > sp->res[1]);
    if(!smooth_status)
    {
        if(ct.verbose && pct.gridpe==0)
            printf("REDUCING   %d   %14.8e  %14.8e\n", sp->istate, sp->res[0],sp->res[1]);
        reduce_it = true;
        //do_mgrid = false;  causes hang with some kpoint distributions
        // This is still a little tricky since if it happens to too many states you can
        // converge to a wrong answer so maybe just leave it off for now so that it won't
        // converge at all.
        //for(int idx = 0;idx <pbasis_noncoll;idx++) tmp_psi_t[idx] = saved_psi[idx];
    }


    /* Now do a multigrid cycle */
    if (do_mgrid )
    {

        for(int is = 0; is < ct.noncoll_factor; is++)
        {

            /* Do multigrid step with solution returned in sg_twovpsi */
            {

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
                        NX_GRID, NY_GRID, NZ_GRID,
                        G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                        dimx, dimy, dimz, ct.boundaryflag);

                MG.mg_prolong (twork_tf, v_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
                CopyAndConvert(sbasis, (mgtype_t *)twork_tf, (convert_type_t *)sg_twovpsi_t);

            }

            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */
            if(ct.use_rmm_diis)
            {
                CPP_pack_stop<CalcType> (sg_twovpsi_t, &work2_t[is*pbasis], dimx, dimy, dimz);
                for(int i=0;i < pbasis_noncoll;i++)
                {
                    rmmres_t[is*pbasis + i] = (saved_psi[is*pbasis+i] -
                                              tmp_psi_t[is*pbasis+i] +
                                              work2_t[is*pbasis+i]);
                }
            }
            else
            {
                // Gradient step
                CPP_pack_stop_axpy<CalcType>(sg_twovpsi_t, &tmp_psi_t[is*pbasis], -fg_step, dimx, dimy, dimz);
            }
        } 

        if(ct.use_rmm_diis)
        {
            sp->dptr->addres(rmmres_t);  // Preconditioned residual

            // If first vcycle use gradient step
            if(vcycle == 0) sp->dptr->lambda = -0.5*fg_step;
            // compute optimized for second vcycle
            if(vcycle >= 2)
            {
                // line minimization to compute rmm-diis lambda
                ApplyHamiltonian<OrbitalType,CalcType> (
                      kptr, sp,
                      sp->istate,
                      rmmres_t,  // preconditioned residual |kr0>
                      hr0_t,     // hamiltonian applied to  |kro>
                      vtot_psi,
                      vxc_psi,
                      nv_t,
                      false);

                // not clear if initial residual r0_t or final residual 
                // res_t works better here
                sp->dptr->compute_lambda(sp->eig[0], ihu_t, r0_t, hr0_t);
            }

            std::vector<OrbitalType> next;
            next = sp->dptr->compute_estimate(); 
            if(next.size() > 0)
            {
                for(int i=0;i < pbasis_noncoll;i++)
                    tmp_psi_t[i] = (CalcType)next[i] + sp->dptr->lambda*rmmres_t[i];
            }
            else
            {
                // Fall back to gradient step
                for(int i=0;i < pbasis_noncoll;i++)
                    tmp_psi_t[i] = saved_psi[i] - fg_step*rmmres_t[i];
            }
        }
    }

    if(potential_acceleration)
        PotentialAcceleration(kptr, sp, vtot_psi, dvtot_psi, tmp_psi_t, saved_psi);

    // Copy single precision orbital back to double precision
    ScatterPsi(G, pbasis_noncoll, sp->istate, tmp_psi_t, kptr->orbital_storage, pct.coalesce_factor);

    p->free(res2_t, pool_blocks);

} // end MgEigState


