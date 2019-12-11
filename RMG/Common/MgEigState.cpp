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


template void MgEigState<double,float>(Kpoint<double> *, State<double> *, double *, double *, double *, double *, int);
template void MgEigState<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, double *, int);
template void MgEigState<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *,
        double *, std::complex<double> *, std::complex<double> *, int);
template void MgEigState<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *,
        double *, std::complex<double> *, std::complex<double> *, int);



template <typename OrbitalType, typename CalcType>
void MgEigState (Kpoint<OrbitalType> *kptr, State<OrbitalType> * sp, double * vtot_psi, double *vxc_psi, OrbitalType *nv, OrbitalType *ns, int vcycle)
{
    RmgTimer RT("Mg_eig");
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
    int eig_post[MAX_MG_LEVELS] = { 0, 0, 4, 4, 4, 4, 4, 4 };

    int potential_acceleration;
    Mgrid MG(L, T);

    int nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;
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
    if ((ct.runflag == RANDOM_START) && (ct.scf_steps < 2)) do_mgrid = false;

    double Zfac = 2.0 * ct.max_zvalence;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);
    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    int sbasis_noncoll = sbasis * ct.noncoll_factor;

    size_t aratio = sizeof(OrbitalType) / sizeof(CalcType);

    // Boost pool only makes a single call to system allocation routines and manages the blocks
    // after that which reduces contention when many threads are running. Automatically frees
    // the allocated memory when it goes out of scope.
    int pool_blocks = 24;
    boost::pool<> p(sbasis_noncoll*aratio*sizeof(CalcType), pool_blocks);
    CalcType *res2_t = (CalcType *)p.ordered_malloc(1);
    CalcType *work2_t = (CalcType *)p.ordered_malloc(4);
    CalcType *work1_t = (CalcType *)p.ordered_malloc(4);
    CalcType *sg_twovpsi_t  =  (CalcType *)p.ordered_malloc(4);
    OrbitalType *saved_psi  = (OrbitalType *)p.ordered_malloc(aratio);
    double *dvtot_psi = (double *)p.ordered_malloc(aratio);
    CalcType *tmp_psi_t  = (CalcType *)p.ordered_malloc(1);
    CalcType *res_t  =  (CalcType *)p.ordered_malloc(1);
    CalcType *twork_t  = (CalcType *)p.ordered_malloc(1);
    OrbitalType *nv_t  = (OrbitalType *)p.ordered_malloc(aratio);

    // Copy double precision psi into correct precison array
    GatherPsi(G, pbasis_noncoll, sp->istate, kptr->orbital_storage, tmp_psi_t);

    // Copy nv into local array
    if(ct.coalesce_states)
        GatherPsi(G, pbasis_noncoll, sp->istate, kptr->nv, nv_t);
    else
        GatherPsi(G, pbasis_noncoll, 0, nv, nv_t);


    // For USPP copy double precision ns into correct precision temp array. For NCPP ns=psi. */
    if(ct.norm_conserving_pp)
    {
        for(int ix=0;ix < pbasis_noncoll;ix++) work1_t[ix] = tmp_psi_t[ix];
    }
    else
    {
        GatherPsi(G, pbasis_noncoll, sp->istate, kptr->ns, work1_t);
    }

    /* Save in res2 */
    for(int ix=0;ix < pbasis_noncoll;ix++) res2_t[ix] = work1_t[ix];

    // Setup some potential acceleration stuff
    potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    if(potential_acceleration) {
        for(int idx = 0;idx <pbasis_noncoll;idx++) saved_psi[idx] = (OrbitalType)tmp_psi_t[idx];
        PotentialAccelerationWait(sp->istate, kptr->nstates, kptr->dvh_skip);
    }

    /* Smoothing cycles */
    for (int cycles = 0; cycles <= nits; cycles++)
    {
        /* Apply Hamiltonian */
        {
            RmgTimer RT1("Mg_eig: apply hamiltonian");
            diag=ApplyHamiltonian<OrbitalType,CalcType> (kptr, tmp_psi_t, work1_t, vtot_psi, vxc_psi, nv);
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
                    eig = 0.3 * eig + 0.7 * sp->oldeig[0];
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
                            ct.eig_parm.sb_step, 2.0*Zfac, 0.0, NULL,
                            NX_GRID, NY_GRID, NZ_GRID,
                            G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                            dimx, dimy, dimz, ct.boundaryflag);

                    MG.mg_prolong (twork_tf, v_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
                    CopyAndConvert(sbasis, (mgtype_t *)twork_tf, (convert_type_t *)sg_twovpsi_t);

                }

                /* The correction is in a smoothing grid so we use this
                 * routine to update the orbital which is stored in a physical grid.
                 */

                t1 = -ct.eig_parm.mg_timestep;
                CPP_pack_stop_axpy<CalcType> (sg_twovpsi_t, &tmp_psi_t[is*pbasis], t1, dimx, dimy, dimz);

            }
            else
            {

                t1 = eig;
                double t5 = diag - Zfac;
                t5 = -1.0 / t5;
                double t4 = ct.eig_parm.gl_step * t5;
                for (int idx = 0; idx <pbasis; idx++)
                {
                    OrbitalType t5 = t4 * (OrbitalType)res_t[idx + is * pbasis];
                    tmp_psi_t[idx + is * pbasis] += t5;
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
        ScatterPsi(G, pbasis_noncoll, sp->istate, tmp_psi_t, kptr->orbital_storage);

} // end MgEigState


