/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#include <complex>
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "BlasWrappers.h"
#include "auxiliary.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "packfuncs.h"
#include "transition.h"


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


extern STATE *states;


static std::mutex vtot_sync_mutex;

template void MgEigState<double,float>(BaseGrid *, TradeImages *, Lattice *, STATE *, int , double *);
template void MgEigState<std::complex<double>, std::complex<float> >(BaseGrid *, TradeImages *, Lattice *, STATE *, int , double *);

template <typename OrbitalType, typename CalcType>
void MgEigState (BaseGrid *G, TradeImages *T, Lattice *L, STATE * sp, int tid, double * vtot_psi)
{

    RmgTimer RT("Mg_eig");

    int idx, cycles, pbasis;
    int nits, sbasis;
    double eig, diag, t1, t2, t3, t4;
    double *work1;
    OrbitalType *nv, *ns;
    double *res;
    double *nvtot_psi;
    int eig_pre[6] = { 0, 3, 6, 2, 2, 2 };
    int eig_post[6] = { 0, 3, 6, 2, 2, 2 };
    int ione = 1;
    int dimx, dimy, dimz, levels, potential_acceleration;
    double hxgrid, hygrid, hzgrid, sb_step;
    Mgrid MG(L, T);

    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;
    dimx = G->get_PX0_GRID(1);
    dimy = G->get_PY0_GRID(1);
    dimz = G->get_PZ0_GRID(1);
    hxgrid = G->get_hxgrid(1);
    hygrid = G->get_hygrid(1);
    hzgrid = G->get_hzgrid(1);
    levels = ct.eig_parm.levels;
    if ((ct.runflag == 0) && (ct.scf_steps < 2)) {
        levels = 0;
    }

    sb_step = 1.0;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;

    /* Grab some memory */
    CalcType *res2_t = new CalcType[2*sbasis];
    CalcType *work2_t = new CalcType[4 * sbasis];
    CalcType *work1_t = new CalcType[4 * sbasis];
    work1 = new double[4 * sbasis];

    CalcType *sg_psi_t = new CalcType[2 * sbasis];
    res = new double[2*sbasis];
    CalcType *sg_twovpsi_t = new CalcType[2 * sbasis];

    CalcType *res2 = new CalcType[2 * sbasis];
    OrbitalType *saved_psi = new OrbitalType[2 * sbasis];
    nvtot_psi = new double[2*sbasis];
    CalcType *tmp_psi_t = new CalcType[2 * sbasis];
    CalcType *res_t = new CalcType[2 * sbasis];

    OrbitalType *tmp_psi = (OrbitalType *)sp->psiR;
    std::complex<double> *kdr = new std::complex<double>[2*sbasis];
    for(int idx = 0;idx < pbasis;idx++) kdr[idx] = 0.0;


    if(ct.eig_parm.mucycles > 1) {
        mix_betaxpsi1(sp);
    }

    /* Get the non-local operator and S acting on psi (nv and ns, respectfully) */
    if(ct.is_gamma) {
        nv = (OrbitalType *)&pct.nv[sp->istate * pbasis]; 
        ns = (OrbitalType *)&pct.ns[sp->istate * pbasis]; 
    }
    else {
        nv = (OrbitalType *)&pct.nv[2*sp->istate * pbasis]; 
        ns = (OrbitalType *)&pct.ns[2*sp->istate * pbasis]; 
    }


    // Copy double precision ns into temp single precision array */
    CopyAndConvert(pbasis, ns, work1_t);

    /*Apply double precision Mehrstellen right hand operator to ns and save in res2 */
    RmgTimer *RT1 = new RmgTimer("Mg_eig: app_cir");
    CPP_app_cir_driver<CalcType> (L, T, work1_t, res2_t, dimx, dimy, dimz, ct.kohn_sham_fd_order);
    delete(RT1);

    // Copy double precision psi into single precison array
    CopyAndConvert(pbasis, tmp_psi, tmp_psi_t);

    // Setup some potential acceleration stuff
    potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    if(potential_acceleration) {
        for(idx = 0;idx <pbasis;idx++) {
            nvtot_psi[idx] = vtot_psi[idx];
            saved_psi[idx] = tmp_psi[idx];
        }
    }




    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        /* Apply Mehrstellen left hand operators */
        RT1 = new RmgTimer("Mg_eig: app_cil");
        diag = CPP_app_cil_driver<CalcType> (L, T, tmp_psi_t, work2_t, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);
        delete(RT1);
        diag = -1.0 / diag;

        // if complex orbitals apply gradient to psi and compute dot products
        if(typeid(OrbitalType) == typeid(std::complex<double>)) {

            CalcType *gx = new CalcType[pbasis];
            CalcType *gy = new CalcType[pbasis];
            CalcType *gz = new CalcType[pbasis];

            CPP_app_grad_driver (L, T, tmp_psi_t, gx, gy, gz, dimx, dimy, dimz, hxgrid, hygrid, hzgrid);

            std::complex<double> I_t(0.0, 1.0);
            for(int idx = 0;idx < pbasis;idx++) {

                kdr[idx] = -I_t * (ct.kp[sp->kidx].kvec[0] * (std::complex<double>)gx[idx] +
                                               ct.kp[sp->kidx].kvec[1] * (std::complex<double>)gy[idx] +
                                               ct.kp[sp->kidx].kvec[2] * (std::complex<double>)gz[idx]);
            }

            delete [] gz;
            delete [] gy;
            delete [] gx;
        }

        // Copy saved application to ns to res
        for(int idx=0;idx < pbasis;idx++) {
            res_t[idx] = res2_t[idx];
        }

        if(potential_acceleration) {
            /* Generate 2 * V * psi */
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, nvtot_psi, (void *)kdr, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);
        }
        else {
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, vtot_psi, (void *)kdr, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);
        }

        /* B operating on 2*V*psi stored in work1 */
        RT1 = new RmgTimer("Mg_eig: app_cir");
        CPP_app_cir_driver<CalcType> (L, T, sg_twovpsi_t, work1_t, dimx, dimy, dimz, ct.kohn_sham_fd_order);
        delete(RT1);

        // Add in non-local which has already had B applied in AppNls
        for(int idx=0;idx < pbasis;idx++) {
            work1_t[idx] += 2.0 * nv[idx];
        }

        for(int idx=0;idx < pbasis;idx++) {
            work1_t[idx] = work1_t[idx] - work2_t[idx];
        }

        /* If this is the first time through compute the eigenvalue */
        if ((cycles == 0) || (potential_acceleration != 0)) 
        {

            ComputeEig(pbasis, tmp_psi_t, work1_t, res_t, &eig);
            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
	     * are not defined yet*/
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

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {


            CalcType f1(TWO * eig);
            for (idx = 0; idx <pbasis; idx++)
            {

                res_t[idx] = f1 * res_t[idx] - work1_t[idx];

            }

            /* Pack the residual data into multigrid array */
            CPP_pack_ptos<CalcType> (sg_psi_t, res_t, dimx, dimy, dimz);

            T->trade_images<CalcType> (sg_psi_t, dimx, dimy, dimz, FULL_TRADE);


            /* Smooth it once and store the smoothed residual in work1 */
            t1 = 145.0;

            CPP_app_smooth<CalcType> (sg_psi_t, work1_t, G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1));


            if(potential_acceleration) {
                t1 = eig - states[0].eig[0];
                t1 = -t1*t1 / 10.0;
            }
            else {
                t1 = 0.0;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            RT1 = new RmgTimer("Mg_eig: mgrid_solv");
            MG.mgrid_solv<CalcType> (sg_twovpsi_t, work1_t, work2_t,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, sb_step, t1,
                        G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                        G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                        G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
            delete(RT1);

            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */

            t1 = -1.0;
            CPP_pack_stop_axpy<CalcType> (sg_twovpsi_t, tmp_psi_t, t1, dimx, dimy, dimz);


        }
        else
        {


            t1 = TWO * eig;
            t2 = ZERO;
            t4 = ct.eig_parm.gl_step * diag;

            for (idx = 0; idx <pbasis; idx++)
            {

                OrbitalType t5;
                t5 = t1 * (OrbitalType)res_t[idx] - (OrbitalType)work1_t[idx];
                //RRRRt2 += t5 * t5;
                tmp_psi_t[idx] += t4 * t5;

            }

            if (cycles == 0)
            {

                t2 = real_sum_all (t2, pct.grid_comm);
                t1 = (double) (ct.psi_nbasis);
                sp->res = ct.hmaxgrid * ct.hmaxgrid * sqrt (t2 / t1) * 0.25;

            }

        }

    }                           /* end for */
#if 0
    if(potential_acceleration) {

        // Save potential used for this orbital and update potential for future orbitals
        for(idx = 0;idx <pbasis;idx++) {
            sp->dvhxc[idx] = nvtot_psi[idx];
        }


        if(ct.potential_acceleration_constant_step > 0.0) {

            t1 = 1.8 * ct.potential_acceleration_constant_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;

            vtot_sync_mutex.lock();
            for(idx = 0;idx <pbasis;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * PI * sp->occupation[0] * tmp_psi_t[idx] * (tmp_psi_t[idx] - (CalcType)saved_psi[idx]);
            }
            vtot_sync_mutex.unlock();

        }

        if(ct.potential_acceleration_poisson_step > 0.0) {

            // construct delta_rho
            for(idx = 0;idx <pbasis;idx++) {
                res_t[idx] = -4.0 * PI * sp->occupation[0] *
                           (tmp_psi_t[idx] - saved_psi[idx]) * (2.0*saved_psi[idx] + (tmp_psi_t[idx] - saved_psi[idx]));
            }

            // zero out solution vector
            for(idx = 0;idx <pbasis;idx++) {
                sg_twovpsi_t[idx] = 0.0;
            }

            /* Pack delta_rho into multigrid array */
            CPP_pack_ptos<CalcType> (sg_psi_t, res_t, dimx, dimy, dimz);
            T->trade_images<CalcType> (sg_psi_t, dimx, dimy, dimz, FULL_TRADE);
            /* Smooth it once and store the smoothed charge in res */
            CPP_app_smooth1<CalcType> (sg_psi_t, res_t, G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1));

            // neutralize cell with a constant background charge
            t2 = 0.0;
            for(idx = 0;idx <pbasis;idx++) {
                t2 += res_t[idx];
            }
            t2 = real_sum_all(t2, pct.grid_comm) / (G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));
            for(idx = 0;idx <pbasis;idx++) {
                res_t[idx] -= t2;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            eig_pre[0] = 2;
            eig_post[0] = 2;
            levels=1;

            RT1 = new RmgTimer("Mg_eig: mgrid_solv");
            MG.mgrid_solv (sg_twovpsi_t, res_t, work2_t,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, G->get_neighbors(), levels, eig_pre, eig_post, 1, 1.0, 0.0,
                        G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                        G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                        G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
            delete(RT1);

            for(idx = 0;idx <pbasis;idx++) {
                res_t[idx] = 0.0;
            }
            CPP_pack_stop_axpy<CalcType> (sg_twovpsi_t, res_t, 1.0, dimx, dimy, dimz);
            t1 = ct.potential_acceleration_poisson_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;
            vtot_sync_mutex.lock();
            for(idx = 0;idx <pbasis;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * res_t[idx];
            }
            vtot_sync_mutex.unlock();
        }


    } // end if
#endif

    // Copy single precision orbital back to double precision
    CopyAndConvert(pbasis, tmp_psi_t, tmp_psi);


    /* Release our memory */
    delete [] kdr;
    delete [] res_t;
    delete [] tmp_psi_t;
    delete [] nvtot_psi;
    delete [] saved_psi;
    delete [] res2;

    delete [] sg_twovpsi_t;
    delete [] res;
    delete [] sg_psi_t;

    delete [] work1;
    delete [] work1_t;
    delete [] work2_t;
    delete [] res2_t;



} // end MgEigState

#include "transition.h"

extern "C" void mg_eig_state(STATE *sp, int tid, double *vtot_psi)
{
    MgEigState<double,float> (Rmg_G, Rmg_T, &Rmg_L, sp, tid, vtot_psi);
}

extern "C" void mg_eig_state_driver (STATE * sp, int tid, double * vtot_psi)
{

        if(ct.is_gamma) {
            MgEigState<double,float> (Rmg_G, Rmg_T, &Rmg_L, sp, tid, vtot_psi);
        }
        else {
            MgEigState<std::complex<double>,std::complex<float> > (Rmg_G, Rmg_T, &Rmg_L, sp, tid, vtot_psi);
        }


}

