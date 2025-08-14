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


template <typename OrbitalType, typename CalcType> void Smoother(
              Kpoint<OrbitalType> *kp, State<OrbitalType> *sp,
              CalcType *u, CalcType *Hu, CalcType *r,
              double *v, double *vxc, double *dinv,
              OrbitalType *nv, CalcType *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin);


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
    double lambda_max=6.0, lambda_min=0.7;   // Non-ideal placeholders for now
    BaseThread *Thread = BaseThread::getBaseThread(0);
    int tid = Thread->get_thread_tid();

    // Save in case needed for variational energy correction term
    sp->feig[0]=sp->eig[0];

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

    double eig=0.0, t1;
    int eig_pre[MAX_MG_LEVELS] = { 0, 8, 8, 8, 8, 8, 8, 8 };
    int eig_post[MAX_MG_LEVELS] = { 0, 8, 8, 8, 8, 8, 8, 8 };
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

    wfobj<double> dinv;
    for(int ix=0;ix < dinv.pbasis;ix++)
        dinv[ix] = 1.0 / std::abs(diag + 2.0*vtot_psi[ix] + 0.5*kptr->kp.kmag);

    //bool is_jacobi = false;
    bool is_jacobi = true;

    Smoother<OrbitalType,CalcType>(kptr, sp,
              tmp_psi_t, work1_t, res_t,
              vtot_psi, vxc_psi, dinv.data(),
              nv_t, res2_t,
              sp->eig[0], 3, is_jacobi, lambda_max, lambda_min);

    /* Now do a multigrid cycle */
    if (do_mgrid )
    {

        for(int is = 0; is < ct.noncoll_factor; is++)
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
    }

    // Post smoothing
    // Once the SCF error gets small increase this to remove
    // interpolation errors from coarse to fine.
    int gl_pst = 0;
    if(ct.scf_accuracy < 1.0e-8) gl_pst=1;
    is_jacobi = true;
    for(int its=0;its < gl_pst;its++)
    {
        Smoother<OrbitalType,CalcType>(kptr, sp,
                  tmp_psi_t, work1_t, res_t,
                  vtot_psi, vxc_psi, dinv.data(),
                  nv_t, res2_t,
                  sp->eig[0], 0, is_jacobi, lambda_max, lambda_min);
    }

    if(potential_acceleration)
        PotentialAcceleration(kptr, sp, vtot_psi, dvtot_psi, tmp_psi_t, saved_psi);

    // Copy single precision orbital back to double precision
    if(freeze_occupied)
        ScatterPsi(G, pbasis_noncoll, sp->istate, tmp_psi_t, kptr->orbital_storage, pct.coalesce_factor);

    p->free(res2_t, pool_blocks);

} // end MgEigState


// Multigrid fine grid smoothing routine. If is_jacobi is true
// it performs order iterations of relaxed Jacobi with omega = 2/3.
// If is_jacobi is false it performs Chebyshev of the requested order.
// If order is zero it does jacobi without a residual calculation on exit.
template <typename OrbitalType, typename CalcType>
void Smoother (Kpoint<OrbitalType> *kptr,
               State<OrbitalType> *sp,
               CalcType *u,        // wavefunction being smoothed
               CalcType *Hu,       // if order>0 on exit holds last Hu
               CalcType *r,        // workspace passed in from caller
               double *v,
               double *vxc,
               double *dinv,
               OrbitalType *nv,
               CalcType *ns,       // only used for USPP
               double &eig, 
               int order,
               bool is_jacobi,
               double lmax,
               double lmin)
{

    const double theta = 0.5*(lmax + lmin);
    const double delta = 0.5*(lmax - lmin);
    const double sigma = theta / delta;

    wfobj<CalcType> p;
    CalcType f1;

    /* Apply Hamiltonian */
    {
        RmgTimer RT1("Mg_eig: apply hamiltonian");
        ApplyHamiltonian<OrbitalType,CalcType> (kptr, sp, sp->istate, u, Hu, v, vxc, nv, false);
    }

    double eig1 = ComputeEig(p.pbasis, u, Hu, ns);
    eig = 0.7*eig1 + 0.3*eig;
    f1 = 2.0*eig;

    for(int i=0;i < p.pbasis;i++)
    {
        //r[i] = f1*ns[i] - 2.0*Hu[i];
        r[i] = f1*u[i] - 2.0*Hu[i];
        p[i] = dinv[i] * r[i];
    }

    double a=0.0, b=0;
    for(int k=0;k < order;k++)
    {
        if (k==0)
        {
            a = 2.0 / theta; b = 0.0;
        }
        else if (k==1)
        {
            a = 2.0 / theta;
            b  = (a * delta * delta) / 4.0;
        }
        else
        {
            const double beta_new = 1.0 / (2.0*sigma - b);
            b = beta_new;
            a = 2.0 * beta_new / delta;
        }

        if(is_jacobi)
        {
            a = 2.0/3.0;
            b = 0.0;
        }

        //if(pct.gridpe==0 && sp->istate==0)printf("QQQQ  %f  %f\n",a, b);

        for(int i=0;i < p.pbasis;i++) u[i] += a*p[i];
        ApplyHamiltonian<OrbitalType,CalcType> (kptr, sp, sp->istate, u, Hu, v, vxc, nv, false);
        eig1 = ComputeEig(p.pbasis, u, Hu, ns);
        eig = 0.7*eig1 + 0.3*eig;
        f1 = 2.0*eig;

        double rsum = 0.0;
        for (int i=0;i<p.pbasis;i++)
        {
            //r[i] = f1*ns[i] - 2.0*Hu[i];
            r[i] = f1*u[i] - 2.0*Hu[i];
            p[i] = dinv[i]*r[i] + b * p[i];
            rsum += std::real(r[i] * std::conj(r[i]));
        }
        GlobalSums (&rsum, 1, pct.coalesced_grid_comm);
        rsum *= get_vel();
    }
}
template void Smoother<double,float>(
              Kpoint<double> *kp, State<double> *sp,
              float *u, float *Hu, float *r,
              double *v, double *vxc, double *dinv,
              double *nv, float *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin);

template void Smoother<double,double>(
              Kpoint<double> *kp, State<double> *sp,
              double *u, double *Hu, double *r,
              double *v, double *vxc, double *dinv,
              double *nv, double *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin);

template void Smoother<std::complex<double>,std::complex<float>>(
              Kpoint<std::complex<double>> *kp, State<std::complex<double>> *sp,
              std::complex<float> *u, std::complex<float> *Hu, std::complex<float> *r,
              double *v, double *vxc, double *dinv,
              std::complex<double> *nv, std::complex<float> *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin);

template void Smoother<std::complex<double>,std::complex<double>>(
              Kpoint<std::complex<double>> *kp, State<std::complex<double>> *sp,
              std::complex<double> *u, std::complex<double> *Hu, std::complex<double> *r,
              double *v, double *vxc, double *dinv,
              std::complex<double> *nv, std::complex<double> *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin);

