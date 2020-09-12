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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include <boost/math/special_functions/erf.hpp>
#include "Atomic.h"

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;


// Generates a set of projectors according to the Bloechl prescription
// using spherical Bessel functions as the basis set.
void InitSemilocalBessel (void)
{

    Atomic A;
    double nlradius = 7.5;
    //ct.max_l = 0;
    int max_nlprojectors = 0;

    for (auto sp = std::begin (Species); sp != std::end (Species); ++sp)
    {

        std::vector<double> work, work2;
        work.resize(sp->rg_points);
        work2.resize(sp->rg_points);
        sp->llbeta.clear();
        sp->nbeta = 0;
        sp->nhtol.clear();
        sp->nhtom.clear();
        sp->indv.clear();
        sp->nh_l2m.clear();
        sp->nhtoj.clear();

        int beta_idx = 0;
        int total_beta = 0;
        for(size_t il = 0; il < sp->dVl_l.size(); il++)
        {
            int lval = sp->dVl_l[il];
            int N = A.CountRoots(lval, nlradius, ct.rhocparm, 4.0*ct.hmingrid);
            total_beta += N;
        }

        double_2d_array ddd0; 
        ddd0.resize(boost::extents[total_beta][total_beta]);
        for(int i = 0;i < total_beta;i++)
            for(int j = 0;j < total_beta;j++)
                ddd0[i][j] = 0.0;

        for(size_t il = 0; il < sp->dVl_l.size(); il++)
        {
            int lval = sp->dVl_l[il];
            int N = A.CountRoots(lval, nlradius, ct.rhocparm, 4.0*ct.hmingrid);
            double_2d_array phi;
            phi.resize(boost::extents[N][sp->rg_points]);

            // Generate the basis set
            double anorm = 1.0 / sqrt(4.0 * PI);
            for(int ib = 0;ib < N;ib++)
            {
                sp->llbeta.emplace_back(lval);
                A.GenBessel(phi.origin() + ib*sp->rg_points, sp->r, nlradius, sp->rg_points, lval, ib);
                for(int idx=0;idx < sp->rg_points;idx++) phi[ib][idx] *= anorm;
            }

            std::vector<double> ci;
            ci.resize(N);
            std::fill(ci.begin(), ci.end(), 0.0);
            // Orthogonalization procedure
            for(int ix=0;ix < N;ix++)
            {   
                for(int idx=0;idx < sp->rg_points;idx++) work[idx] = phi[ix][idx];
                for(int jx=0;jx < ix;jx++)
                {
                    for(int idx=0;idx < sp->rg_points;idx++) work2[idx] = phi[jx][idx]*sp->dVl[il][idx]*phi[ix][idx];
                    double t1 = radint1 (work2.data(), sp->r, sp->rab, sp->rg_points);
                    for(int idx=0;idx < sp->rg_points;idx++) work[idx] -= phi[jx][idx] * ci[jx] * t1;
                }
                for(int idx=0;idx < sp->rg_points;idx++) phi[ix][idx] = work[idx];
                for(int idx=0;idx < sp->rg_points;idx++) work[idx] = phi[ix][idx]*sp->dVl[il][idx]*phi[ix][idx];
                ci[ix] = radint1 (work.data(), sp->r, sp->rab, sp->rg_points);
                ci[ix] = 1.0 / ci[ix];
            }

            // Generate new beta
            for(int ix=0;ix < N;ix++)
            {
                sp->beta.emplace_back(new double[sp->rg_points]);
                for(int idx=0;idx < sp->rg_points;idx++) sp->beta[sp->nbeta][idx] = ci[ix] * sp->dVl[il][idx] * phi[ix][idx];
                sp->nbeta++;
            }

            for(int i=0;i < N;i++) ddd0[i+beta_idx][i+beta_idx] = 1.0 / ci[i];
            beta_idx += N;

        }

        // Count up projectors in order to resize arrays correctly
        int ihmax = 0;
        for (int j = 0; j < sp->nbeta; j++)
        {   
            int l = sp->llbeta[j];
            for (int k = 0; k < 2 * l + 1; k++)
            {   
                ++ihmax;
            }
        }

        for (int j = 0; j < ihmax; j++)
        {   
            sp->nhtol.push_back(0);
            sp->nhtom.push_back(0);
            sp->indv.push_back(0);
            sp->nh_l2m.push_back(0);
            sp->nhtoj.push_back(0.0);
        }
      
        int ih = 0;
        for (int j = 0; j < sp->nbeta; j++)
        {
            int l = sp->llbeta[j];
            for (int k = 0; k < 2 * l + 1; k++)
            {
                sp->nhtol[ih] = l;
                sp->nhtom[ih] = k;
                sp->indv[ih] = j;
                sp->nh_l2m[ih] = l*l + k;
                if(sp->is_spinorb) sp->nhtoj[ih] = sp->jjbeta[j];
                ++ih;
            }
        }
        sp->nh = ih;
        if (ih > max_nlprojectors)
            max_nlprojectors = ih;

     
        sp->ddd0.resize(boost::extents[ih][ih]);
        sp->qqq.resize(boost::extents[ih][ih]);
        for (int j = 0; j < ih; j++)
        {
            for (int k = 0; k < ih; k++)
            {
                sp->ddd0[j][k] = ddd0[sp->indv[j]][sp->indv[k]];
//                sp->qqq[j][k] = qqq[sp->indv[j]][sp->indv[k]];
                /*                          sp->ddd[j][k]=ddd[indv[j]][indv[k]]; */
                if(sp->is_spinorb) continue;
                if ((sp->nhtol[j] != sp->nhtol[k]) || (sp->nhtom[j] != sp->nhtom[k]))
                {
                    sp->ddd0[j][k] = 0.0;
                    sp->qqq[j][k] = 0.0;
                    /*                          sp->ddd[j][k]=0.0; */
                }
            }
        }

        // Set the maximum number of non-local projecters needed
        ct.max_nl = std::max(ct.max_nl, max_nlprojectors);

    }
}
