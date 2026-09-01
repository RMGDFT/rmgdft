/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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


template <typename OrbitalType, typename CalcType> void mgsmoother(
              Kpoint<OrbitalType> *kp, State<OrbitalType> *sp,
              CalcType *u, CalcType *Hu, CalcType *r, 
              CalcType *iHu, CalcType *ir,
              double *v, double *vxc, double *dinv,
              OrbitalType *nv, CalcType *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin, int vcycle);


// Multigrid fine grid smoothing routine. If is_jacobi is true
// it performs order iterations of relaxed Jacobi with omega = 2/3.
// If is_jacobi is false it performs Chebyshev of the requested order.
// If order is zero it does jacobi without a residual calculation on exit.
template <typename OrbitalType, typename CalcType>
void mgsmoother (Kpoint<OrbitalType> *kptr,
               State<OrbitalType> *sp,
               CalcType *u,        // wavefunction being smoothed
               CalcType *Hu,       // if order>0 on exit holds last Hu
               CalcType *r,        // on exit holds last residual
               CalcType *iHu,      // on exit holds first Hu
               CalcType *ir,       // on exit holds first residual
               double *v,
               double *vxc,
               double *dinv,
               OrbitalType *nv,
               CalcType *ns,       // only used for USPP
               double &eig, 
               int order,
               bool is_jacobi,
               double lmax,
               double lmin,
               int vcycle)
{
    const double theta = 0.5*(lmax + lmin);
    const double delta = 0.5*(lmax - lmin);
    const double sigma = theta / delta;
    wfobj<CalcType> pp;
    std::vector<CalcType> p;
    int pbasis = pp.pbasis * ct.noncoll_factor * pct.coalesce_factor;
    p.resize(pbasis);

    ApplyHamiltonian<OrbitalType,CalcType> (kptr, sp, sp->istate, u, Hu, v, vxc, nv, false);
    for (int i=0;i<pbasis;i++) iHu[i] = Hu[i];

    double eig1 = ComputeEig(pbasis, u, Hu, ns);
    auto mixeig = [](double &eig, double &eig1) {
        if((ct.scf_steps && (ct.scf_accuracy < 1.0e-4)))
            eig = eig1;
        else
            eig = 0.7*eig1 + 0.3*eig;
    };
    mixeig(eig, eig1);
    CalcType f1 = 2.0*eig;

    double rsum[2] = {0.0, 0.0};
    for(int i=0;i < pbasis;i++)
    {
        if(ct.norm_conserving_pp)
            r[i] = f1*u[i] - 2.0*Hu[i];
        else
            r[i] = f1*ns[i] - 2.0*Hu[i];
        ir[i] = r[i];
        p[i] = dinv[i] * r[i];
        rsum[0] += std::norm(r[i]);
        rsum[1] += std::norm(u[i]);
    }
    GlobalSums (rsum, 2, pct.coalesced_grid_comm);
    sp->res[0] = rsum[0]*get_vel();
    double norm = rsum[1]*get_vel();
    norm = 1.0 / sqrt(norm);
    if(ct.norm_conserving_pp)
        for(int i=0;i < pbasis;i++) u[i] *= norm;

    //if(pct.gridpe==0 && pct.spinpe==0)
    //    printf("ZZZZ  %d  %14.8e  %14.8e\n",sp->istate, rsum[1]*get_vel(), sp->res[0]);

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

        rsum[0] = 0.0;
        for(int i=0;i < pbasis;i++)
        {
            u[i] += a*p[i];
            rsum[0] += std::norm(u[i]);
        }
        GlobalSums (rsum, 1, pct.coalesced_grid_comm);
        norm = rsum[0]*get_vel();
        norm = 1.0 / sqrt(norm);
        if(ct.norm_conserving_pp)
            for(int i=0;i < pbasis;i++) u[i] *= norm;

        ApplyHamiltonian<OrbitalType,CalcType> (kptr, sp, sp->istate, u, Hu, v, vxc, nv, false);
        eig1 = ComputeEig(pbasis, u, Hu, ns);
        mixeig(eig, eig1);
        f1 = 2.0*eig;

        rsum[0] = 0.0, rsum[1] = 0.0;
        for (int i=0;i<pbasis;i++)
        {
            if(ct.norm_conserving_pp)
                r[i] = f1*u[i] - 2.0*Hu[i];
            else
                r[i] = f1*ns[i] - 2.0*Hu[i];
            p[i] = dinv[i]*r[i] + b * p[i];
            rsum[0] += std::norm(r[i]);
        }

        GlobalSums (rsum, 1, pct.coalesced_grid_comm);
        sp->res[k+1] = rsum[0]*get_vel();
        //if(order >0 && pct.gridpe==0)printf("ZZZZ  %d  %14.8e  \n",sp->istate,rsum[1]*get_vel(), sp->res[k+1]);
    }
}

template void mgsmoother<double,float>(
              Kpoint<double> *kp, State<double> *sp,
              float *u, float *Hu, float *r, float *iHu, float *ir,
              double *v, double *vxc, double *dinv,
              double *nv, float *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin, int vcycle);

template void mgsmoother<double,double>(
              Kpoint<double> *kp, State<double> *sp,
              double *u, double *Hu, double *r, double *iHu, double *ir,
              double *v, double *vxc, double *dinv,
              double *nv, double *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin, int vcycle);

template void mgsmoother<std::complex<double>,std::complex<float>>(
              Kpoint<std::complex<double>> *kp, State<std::complex<double>> *sp,
              std::complex<float> *u, std::complex<float> *Hu,
              std::complex<float> *r, std::complex<float> *iHu, std::complex<float> *ir,
              double *v, double *vxc, double *dinv,
              std::complex<double> *nv, std::complex<float> *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin, int vcycle);

template void mgsmoother<std::complex<double>,std::complex<double>>(
              Kpoint<std::complex<double>> *kp, State<std::complex<double>> *sp,
              std::complex<double> *u, std::complex<double> *Hu,
              std::complex<double> *r, std::complex<double> *iHu, std::complex<double> *ir,
              double *v, double *vxc, double *dinv,
              std::complex<double> *nv, std::complex<double> *ns,
              double &eig, int order, bool is_jacobi, double lmax, double lmin, int vcycle);

