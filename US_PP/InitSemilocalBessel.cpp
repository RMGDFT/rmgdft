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

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;


// Generates a set of radial projectors according to the Bloechl prescription
// using spherical Bessel functions as the basis set. These are later passed
// to the standard radial setup methods.
void SPECIES::InitSemilocalBessel (void)
{

    double nlradius = 3.0;
    int nproj=10;
    double mingrid = nlradius / nproj;
    int max_nlprojectors = 0;

    std::vector<double> work(this->rg_points, 0.0), work2(this->rg_points, 0.0);
    this->llbeta.clear();
    this->beta.clear();
    this->nhtol.clear();
    this->nhtom.clear();
    this->indv.clear();
    this->nh_l2m.clear();
    this->nhtoj.clear();
    this->nbeta = 0;
    int total_beta = 0;

    for(size_t il = 0; il < this->dVl_l.size(); il++)
    {
        int lval = this->dVl_l[il];
        int N = CountRoots(lval, nlradius, ct.rhocparm, mingrid);
        total_beta += 2*N;
    }

    // ddd0 will hold the normalization coefficients for the projectors
    double_2d_array ddd0; 
    ddd0.resize(boost::extents[total_beta][total_beta]);
    std::fill(ddd0.origin(), ddd0.origin() + total_beta*total_beta, 0.0);

    for(size_t il = 0; il < this->dVl_l.size(); il++)
    {
        int lval = this->dVl_l[il];
        int N = CountRoots(lval, nlradius, ct.rhocparm, mingrid);
        double_2d_array phi;
        phi.resize(boost::extents[N][this->rg_points]);

        // Generate the basis set
        double anorm = sqrt(1.0 / 4.0 * PI);
        for(int ib = 0;ib < N;ib++)
        {
            this->llbeta.emplace_back(lval);
            GenBessel(phi.origin() + ib*this->rg_points, this->r, nlradius, this->rg_points, lval, ib);
            for(int idx=0;idx < this->rg_points;idx++) phi[ib][idx] *= anorm;
            for(int idx=0;idx < this->rg_points;idx++) if(this->r[idx] > nlradius) phi[ib][idx] *= 0.0;
        }

        std::vector<double> ci(N, 0.0);

        // Orthogonalization procedure
        for(int ix=0;ix < N;ix++)
        {   
            for(int idx=0;idx < this->rg_points;idx++) work[idx] = phi[ix][idx];
            for(int jx=0;jx < ix;jx++)
            {
                for(int idx=0;idx < this->rg_points;idx++) work2[idx] = phi[jx][idx]*this->dVl[il][idx]*phi[ix][idx];
                double t1 = radint1 (work2.data(), this->r, this->rab, this->rg_points);
                for(int idx=0;idx < this->rg_points;idx++) work[idx] -= phi[jx][idx] * ci[jx] * t1;
            }
            for(int idx=0;idx < this->rg_points;idx++) phi[ix][idx] = work[idx];
            for(int idx=0;idx < this->rg_points;idx++) work[idx] = phi[ix][idx]*this->dVl[il][idx]*phi[ix][idx];
            ci[ix] = radint1 (work.data(), this->r, this->rab, this->rg_points);
            ci[ix] = 1.0 / ci[ix];
        }
        // Generate new beta
        for(int ix=0;ix < N;ix++)
        {
            this->beta.emplace_back(new double[this->rg_points]);
            for(int idx=0;idx < this->rg_points;idx++) this->beta[this->nbeta][idx] = ci[ix] * this->dVl[il][idx] * phi[ix][idx];
            ddd0[this->nbeta][this->nbeta] = 1.0 / ci[ix];
            this->nbeta++;
        }
    }

    // Count up projectors in order to resize arrays correctly
    int ihmax = 0;
    for (int j = 0; j < this->nbeta; j++)
    {   
        int l = this->llbeta[j];
        for (int k = 0; k < 2 * l + 1; k++)
        {   
            ++ihmax;
        }
    }

    for (int j = 0; j < ihmax; j++)
    {   
        this->nhtol.push_back(0);
        this->nhtom.push_back(0);
        this->indv.push_back(0);
        this->nh_l2m.push_back(0);
        this->nhtoj.push_back(0.0);
    }
  
    int ih = 0;
    for (int j = 0; j < this->nbeta; j++)
    {
        int l = this->llbeta[j];
        for (int k = 0; k < 2 * l + 1; k++)
        {
            this->nhtol[ih] = l;
            this->nhtom[ih] = k;
            this->indv[ih] = j;
            this->nh_l2m[ih] = l*l + k;
            if(this->is_spinorb) this->nhtoj[ih] = this->jjbeta[j];
            ++ih;
        }
    }
    this->nh = ih;
    if (ih > max_nlprojectors)
        max_nlprojectors = ih;

 
    this->ddd0.resize(boost::extents[ih][ih]);
    this->qqq.resize(boost::extents[ih][ih]);
    for (int j = 0; j < ih; j++)
    {
        for (int k = 0; k < ih; k++)
        {
            this->ddd0[j][k] = ddd0[this->indv[j]][this->indv[k]];
//                this->qqq[j][k] = qqq[this->indv[j]][this->indv[k]];
            /*                          this->ddd[j][k]=ddd[indv[j]][indv[k]]; */
            if(this->is_spinorb) continue;
            if ((this->nhtol[j] != this->nhtol[k]) || (this->nhtom[j] != this->nhtom[k]))
            {
                this->ddd0[j][k] = 0.0;
                this->qqq[j][k] = 0.0;
                /*                          this->ddd[j][k]=0.0; */
            }
        }
    }

    // Set the maximum number of non-local projecters needed
    ct.max_nl = std::max(ct.max_nl, max_nlprojectors);

}
