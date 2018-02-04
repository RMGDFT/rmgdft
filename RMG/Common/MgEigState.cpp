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
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "packfuncs.h"
#include "transition.h"
#include "RmgParallelFft.h"
#include "BaseThread.h"


void CopyAndConvert(int n, double *A, float *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (float)A[idx];
}

void CopyAndConvert(int n, double *A, double *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = A[idx];
}

void CopyAndConvert(int n, std::complex<double> *A, std::complex<float> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<float>)A[idx];
}
void CopyAndConvert(int n, std::complex<double> *A, std::complex<double> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<double>)A[idx];
}

void CopyAndConvert(int n, float *A, double *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (double)A[idx];
}

void CopyAndConvert(int n, std::complex<float> *A, std::complex<double> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<double>)A[idx];
}

void ComputeEig(int n, float *A, float *B, float *D, double *rval)
{
    double s1[2];
    s1[0] = 0.0;
    s1[1] = 0.0;

    for(int idx = 0;idx < n;idx++) {
        s1[0] += (double)A[idx] * (double)B[idx];
        s1[1] += (double)A[idx] * (double)D[idx];
    }

    int length = 2;
    GlobalSums (s1, length, pct.grid_comm);

    *rval = (s1[0] / (2.0 * s1[1]));

}

void ComputeEig(int n, double *A, double *B, double *D, double *rval)
{
    double s1[2];
    s1[0] = 0.0;
    s1[1] = 0.0;

    for(int idx = 0;idx < n;idx++) {
        s1[0] += (double)A[idx] * (double)B[idx];
        s1[1] += (double)A[idx] * (double)D[idx];
    }

    int length = 2;
    GlobalSums (s1, length, pct.grid_comm);

    *rval = (s1[0] / (2.0 * s1[1]));

}

void ComputeEig(int n, std::complex<float> *A, std::complex<float> *B, std::complex<float> *D, double *rval)
{
    double s1[4];
    s1[0] = 0.0;
    s1[1] = 0.0;
    s1[2] = 0.0;
    s1[3] = 0.0;

    for(int idx = 0;idx < n;idx++) {
        s1[0] = s1[0] + ((double)std::real(A[idx]) * (double)std::real(B[idx]));
        s1[1] = s1[1] + ((double)std::imag(A[idx]) * (double)std::imag(B[idx]));
        s1[2] = s1[2] + ((double)std::real(A[idx]) * (double)std::real(D[idx]));
        s1[3] = s1[3] + ((double)std::imag(A[idx]) * (double)std::imag(D[idx]));
    }

    int length = 4;
    GlobalSums (s1, length, pct.grid_comm);
    *rval = ((s1[0] + s1[1]) / (2.0 * (s1[2] + s1[3])));

}
void ComputeEig(int n, std::complex<double> *A, std::complex<double> *B, std::complex<double> *D, double *rval)
{
    double s1[4];
    s1[0] = 0.0;
    s1[1] = 0.0;
    s1[2] = 0.0;
    s1[3] = 0.0;

    for(int idx = 0;idx < n;idx++) {
        s1[0] = s1[0] + ((double)std::real(A[idx]) * (double)std::real(B[idx]));
        s1[1] = s1[1] + ((double)std::imag(A[idx]) * (double)std::imag(B[idx]));
        s1[2] = s1[2] + ((double)std::real(A[idx]) * (double)std::real(D[idx]));
        s1[3] = s1[3] + ((double)std::imag(A[idx]) * (double)std::imag(D[idx]));
    }

    int length = 4;
    GlobalSums (s1, length, pct.grid_comm);
    *rval = ((s1[0] + s1[1]) / (2.0 * (s1[2] + s1[3])));

}

extern std::mutex vtot_sync_mutex;

template void MgEigState<double,float>(Kpoint<double> *, State<double> *, double *, double *, double *, int);
template void MgEigState<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, int);
template void MgEigState<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *, std::complex<double> *, std::complex<double> *, int);
template void MgEigState<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *, std::complex<double> *, std::complex<double> *, int);



template <typename OrbitalType, typename CalcType>
void MgEigState (Kpoint<OrbitalType> *kptr, State<OrbitalType> * sp, double * vtot_psi, OrbitalType *nv, OrbitalType *ns, int vcycle)
{
    RmgTimer RT("Mg_eig");
    bool freeze_occupied = true;

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

    double eig=0.0, diag, t1, t2, t4;
    int eig_pre[MAX_MG_LEVELS] = { 0, 8, 8, 20, 20, 20, 20, 20 };
    int eig_post[MAX_MG_LEVELS] = { 0, 2, 2, 2, 2, 2, 2, 2 };

    int potential_acceleration;
    Mgrid MG(L, T);

    int nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;
    int dimx = G->get_PX0_GRID(1) * pct.coalesce_factors[0][0];
    int dimy = G->get_PY0_GRID(1) * pct.coalesce_factors[0][1];
    int dimz = G->get_PZ0_GRID(1) * pct.coalesce_factors[0][2];
    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);
    int levels = ct.eig_parm.levels;
    bool do_mgrid = true;
    if ((ct.runflag == RANDOM_START) && (ct.scf_steps < 2)) do_mgrid = false;

    double Zfac = 2.0 * ct.max_zvalence;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    size_t aratio = sizeof(OrbitalType) / sizeof(CalcType);

    // Boost pool only makes a single call to system allocation routines and manages the blocks
    // after that which reduces contention when many threads are running. Automatically frees
    // the allocated memory when it goes out of scope.
    boost::pool<> p(sbasis*aratio*sizeof(CalcType), 18);
    CalcType *res2_t = (CalcType *)p.ordered_malloc(1);
    CalcType *work2_t = (CalcType *)p.ordered_malloc(2);
    CalcType *work1_t = (CalcType *)p.ordered_malloc(2);
    CalcType *sg_twovpsi_t  =  (CalcType *)p.ordered_malloc(2);
    OrbitalType *saved_psi  = (OrbitalType *)p.ordered_malloc(aratio);
    double *nvtot_psi = (double *)p.ordered_malloc(aratio);
    CalcType *tmp_psi_t  = (CalcType *)p.ordered_malloc(1);
    CalcType *res_t  =  (CalcType *)p.ordered_malloc(1);
    CalcType *twork_t  = (CalcType *)p.ordered_malloc(1);

    OrbitalType *tmp_psi = (OrbitalType *)sp->psi;
    std::complex<double> *kdr = NULL;
    if(typeid(OrbitalType) == typeid(std::complex<double>)) kdr = new std::complex<double>[2*sbasis]();

    // Copy double precision ns into temp single precision array */
    CopyAndConvert(pbasis, ns, work1_t);

    /*Apply double precision Mehrstellen right hand operator to ns and save in res2 */
    {
        RmgTimer RT1("Mg_eig: apply B operator");
        ApplyBOperator<CalcType> (work1_t, res2_t, "Coarse");
    }

    // Copy double precision psi into single precison array
    CopyAndConvert(pbasis, tmp_psi, tmp_psi_t);

    // Setup some potential acceleration stuff
    potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    if(potential_acceleration) {
        vtot_sync_mutex.lock();
        for(int idx = 0;idx <pbasis;idx++) nvtot_psi[idx] = vtot_psi[idx];
        vtot_sync_mutex.unlock();
        for(int idx = 0;idx <pbasis;idx++) saved_psi[idx] = tmp_psi[idx];
    }

    /* Smoothing cycles */
    for (int cycles = 0; cycles <= nits; cycles++)
    {

        /* Apply Mehrstellen left hand operators */
        {
            RmgTimer RT1("Mg_eig: apply A operator");
            diag = ApplyAOperator<CalcType> (tmp_psi_t, work2_t, "Coarse");
        }

        // if complex orbitals apply gradient to psi and compute dot products
        {
            RmgTimer RT1("Mg_eig: apply grad");
            if(typeid(OrbitalType) == typeid(std::complex<double>)) {

                CalcType *gx = new CalcType[pbasis];
                CalcType *gy = new CalcType[pbasis];
                CalcType *gz = new CalcType[pbasis];

                ApplyGradient (tmp_psi_t, gx, gy, gz, APP_CI_EIGHT, "Coarse");

                std::complex<double> I_t(0.0, 1.0);
                for(int idx = 0;idx < pbasis;idx++) {

                    kdr[idx] = -I_t * (kptr->kvec[0] * (std::complex<double>)gx[idx] +
                                                   kptr->kvec[1] * (std::complex<double>)gy[idx] +
                                                   kptr->kvec[2] * (std::complex<double>)gz[idx]);
                }

                delete [] gz;
                delete [] gy;
                delete [] gx;
            }
        }

        // Copy saved application to ns to res
        for(int idx=0;idx < pbasis;idx++) res_t[idx] = res2_t[idx];

        if(potential_acceleration) {
            /* Generate 2 * V * psi */
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, nvtot_psi, (void *)kdr, kptr->kmag, dimx, dimy, dimz);
        }
        else {
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, vtot_psi, (void *)kdr, kptr->kmag, dimx, dimy, dimz);
        }

        /* B operating on 2*V*psi stored in work1 */
        {
            RmgTimer RT1("Mg_eig: apply B operator");
            ApplyBOperator<CalcType> (sg_twovpsi_t, work1_t, "Coarse");
        }

        // Add in non-local which has already had B applied in AppNls
        for(int idx=0;idx < pbasis;idx++) work1_t[idx] += 2.0 * nv[idx];

        for(int idx=0;idx < pbasis;idx++) {
            work1_t[idx] = work1_t[idx] - work2_t[idx];
        }


        /* If this is the first time through compute the eigenvalue */
        if ((cycles == 0) || (potential_acceleration != 0) || (using_davidson && (cycles == 0))) 
        {
            ComputeEig(pbasis, tmp_psi_t, work1_t, res_t, &eig);
            // Save this for variational energy correction
            if((cycles == 0) && (vcycle == 0)) sp->feig[0]=eig;


            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
	     * are not defined yet*/
            if(freeze_occupied) {

                if ((ct.diag == 1) && (potential_acceleration == 0) && (ct.scf_steps < ct.end_diag))
                {
                    if (ct.scf_steps == 0)
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
                    if(ct.scf_steps == 0) {
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
        CalcType f1(TWO * eig);
        for (int idx = 0; idx <pbasis; idx++) res_t[idx] = f1 * res_t[idx] - work1_t[idx];


        /* Now either smooth the wavefunction or do a multigrid cycle */
        if ((cycles == ct.eig_parm.gl_pre) && do_mgrid)
        {

            /* Pack the residual data into multigrid array */
            CPP_pack_ptos<CalcType> (work1_t, res_t, dimx, dimy, dimz);
            T->trade_images (work1_t, dimx, dimy, dimz, FULL_TRADE);


            /* Do multigrid step with solution returned in sg_twovpsi */
            {
                RmgTimer RT1("Mg_eig: mgrid_solv");
                int ixoff, iyoff, izoff;
                int dx2 = MG.MG_SIZE (dimx, 0, G->get_NX_GRID(1), G->get_PX_OFFSET(1), G->get_PX0_GRID(1), &ixoff, ct.boundaryflag);
                int dy2 = MG.MG_SIZE (dimy, 0, G->get_NY_GRID(1), G->get_PY_OFFSET(1), G->get_PY0_GRID(1), &iyoff, ct.boundaryflag);
                int dz2 = MG.MG_SIZE (dimz, 0, G->get_NZ_GRID(1), G->get_PZ_OFFSET(1), G->get_PZ0_GRID(1), &izoff, ct.boundaryflag);
    
		if((dx2 < 0) || (dy2 < 0) || (dz2 < 0)) {
		    printf("Multigrid error: Grid cannot be coarsened. Most likely the current grid is not divisable by 2 or 4. It is recommended to use grid that is, at minimum, divisable by 4. The current grid is %d %d %d" , G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1));
                    exit(0);
		}


                // We use a residual correction multigrid scheme where the right hand side is the residual
                // so single precision is adequate for the correction since the errors from lower precision
                // will be approximately 7 decimal digits smaller than the original error we athe errors from lower precision
                // will be approximately 7 decimal digits smaller than the original error we are correcting for
                if(typeid(CalcType) == typeid(double))
                {
                    float *v_mat = (float *)&sg_twovpsi_t[sbasis];
                    float *f_mat = (float *)&work1_t[sbasis];
                    float *twork_tf = (float *)twork_t;
                    for(int idx = 0;idx < sbasis;idx++) twork_tf[idx] = std::real(work1_t[idx]);
                    MG.mg_restrict<float> (twork_tf, f_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

                    MG.mgrid_solv<float> (v_mat, f_mat, (float *)work2_t,
                                dx2, dy2, dz2, 2.0*hxgrid, 2.0*hygrid, 2.0*hzgrid, 
                                1, G->get_neighbors(), levels, eig_pre, eig_post, 1, 
                                ct.eig_parm.sb_step, 2.0*Zfac, 0.0, NULL,
                                G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                                G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                                G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
                
                    MG.mg_prolong<float> (twork_tf, v_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
                    for(int idx = 0;idx < sbasis;idx++) sg_twovpsi_t[idx] = std::real(twork_tf[idx]);

                }
                else
                {
                    CalcType *v_mat = &sg_twovpsi_t[sbasis];
                    CalcType *f_mat = &work1_t[sbasis];
                    MG.mg_restrict (work1_t, f_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

                    MG.mgrid_solv<CalcType> (v_mat, f_mat, work2_t,
                                dx2, dy2, dz2, 2.0*hxgrid, 2.0*hygrid, 2.0*hzgrid, 
                                1, G->get_neighbors(), levels, eig_pre, eig_post, 1, 
                                ct.eig_parm.sb_step, 2.0*Zfac, 0.0, NULL,
                                G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                                G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                                G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
                
                    MG.mg_prolong (sg_twovpsi_t, v_mat, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
                }
            }

            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */

            t1 = -ct.eig_parm.mg_timestep;
            CPP_pack_stop_axpy<CalcType> (sg_twovpsi_t, tmp_psi_t, t1, dimx, dimy, dimz);

        }
        else
        {

            t1 = TWO * eig;
            t2 = ZERO;
            double t5 = diag - Zfac;
            t5 = -1.0 / t5;
            t4 = ct.eig_parm.gl_step * t5;
            for (int idx = 0; idx <pbasis; idx++)
            {
                t2 += std::norm(res_t[idx]);
                OrbitalType t5 = t4 * (OrbitalType)res_t[idx];
                tmp_psi_t[idx] += t5;
            }

            if (cycles == 0)
            {

                // If occupied orbitals are frozen we compute residuals 
                if(freeze_occupied) {
                    //GlobalSums (&t2, 1, pct.grid_comm);
                    //t2 = RmgSumAll (t2, pct.grid_comm);
                    //t1 = (double) (ct.psi_nbasis);
                    //sp->res = sqrt (t2 / t1);
//                    if(pct.imgpe == 0) std::cout << "Orbital " << sp->istate << " residual = " << sp->res << std::endl;
                }

            }

        }

    }                           /* end for */

    if(potential_acceleration)
        PotentialAcceleration(kptr, sp, vtot_psi, nvtot_psi, tmp_psi_t, saved_psi);

    // Copy single precision orbital back to double precision
    if(freeze_occupied)
        CopyAndConvert(pbasis, tmp_psi_t, tmp_psi);

    /* Release our memory */
    if(typeid(OrbitalType) == typeid(std::complex<double>)) delete [] kdr;

} // end MgEigState

