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
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include <complex>
#include "../Headers/prototypes.h"
#include "RmgParallelFft.h"
#include "GlobalSums.h"


template void GetAugRho<double>(Kpoint<double> **, double *);
template void GetAugRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template <typename KpointType> void GetAugRho(Kpoint<KpointType> **Kpts, double *augrho)
{

    int max_product = (ct.max_nl + 1) * ct.max_nl / 2;
    double *qtpr;

    int pbasis = Kpts[0]->G->get_P0_BASIS(Kpts[0]->G->default_FG_RATIO);
    for(int idx = 0;idx < pbasis;idx++)
        augrho[idx] = 0.0;


    if(!ct.norm_conserving_pp) {

        double *product = new double[max_product];
        KpointType *sint = new KpointType[2 * ct.max_nl];

        for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
        {
            int gion = pct.nonloc_ions_list[ion];
            
            if (pct.Qidxptrlen[gion])
            {
                
                ION *iptr = &ct.ions[gion];
           
                int nh = ct.sp[iptr->species].nh;
                
                int *ivec = pct.Qindex[gion];
                int ncount = pct.Qidxptrlen[gion];
                double *qnmI = pct.augfunc[gion];

                for (int i=0; i < max_product; i++)
                    product[i] = 0.0;

                for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
                {

                    //STATE *sp = ct.kp[kpt].kstate;
                    /* Loop over states and accumulate charge */
                    for (int istate = 0; istate < ct.num_states; istate++)
                    {
                        double t1 = Kpts[kpt]->Kstates[istate].occupation[0] * Kpts[kpt]->kweight;

                        for (int i = 0; i < ct.max_nl; i++)
                        {
                            sint[i] = Kpts[kpt]->newsint_local[istate*pct.num_nonloc_ions*ct.max_nl + ion * ct.max_nl + i];
                        }               /*end for i */

                        int idx = 0;
                        for (int i = 0; i < nh; i++)
                        {
                            for (int j = i; j < nh; j++)
                            {

                                if(i == j) {

                                        product[idx] += t1 * (std::real(sint[i]) * std::real(sint[j]) + std::imag(sint[i]) * std::imag(sint[j]));

                                }
                                else {

                                        product[idx] += 2.0 * t1 * (std::real(sint[i]) * std::real(sint[j]) + std::imag(sint[i]) * std::imag(sint[j]));

                                }
                                idx++;
                            }           /*end for j */
                        }               /*end for i */
                    }                   /*end for istate */
                }                       /*end for kpt */


                int idx = 0;
                for (int i = 0; i < nh; i++)
                {
                    for (int j = i; j < nh; j++)
                    {
                        qtpr = qnmI + idx * ncount;
                        for (int icount = 0; icount < ncount; icount++)
                        {
                            augrho[ivec[icount]] += qtpr[icount] * product[idx];
                        }           /*end for icount */
                        idx++;
                    }               /*end for j */
                }                   /*end for i */


            }                       /*end if */

        }                           /*end for ion */

        GlobalSums(augrho, pbasis, pct.kpsub_comm);
        symmetrize_rho (augrho);

        delete [] sint;
        delete [] product;

    }

}
