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
#include "RmgParallelFft.h"
#include "GlobalSums.h"


void betapsi_so(std::complex<double> *sint, std::complex<double> *sint_tem, SPECIES &sp);

template void GetAugRho<double>(Kpoint<double> **, double *);
template void GetAugRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template <typename KpointType> void GetAugRho(Kpoint<KpointType> **Kpts, double *augrho)
{

    int max_product = (ct.max_nl + 1) * ct.max_nl / 2;

    int pbasis = Kpts[0]->G->get_P0_BASIS(Kpts[0]->G->default_FG_RATIO);
    int num_nonloc_ions = Kpts[0]->BetaProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = Kpts[0]->BetaProjector->get_nonloc_ions_list();

    int factor = ct.noncoll_factor * ct.noncoll_factor;
    for(int idx = 0;idx < pbasis*factor;idx++) augrho[idx] = 0.0;


    if(ct.norm_conserving_pp) return;

    double *product = new double[max_product * factor];
    KpointType *sint = new KpointType[2 * ct.max_nl];

    std::complex<double> *sint_tem = new std::complex<double>[2*ct.max_nl];

    for (int ion = 0; ion < num_nonloc_ions; ion++)
    {
        int gion = nonloc_ions_list[ion];

        if (Atoms[gion].Qindex.size())
        {

            ION *iptr = &Atoms[gion];

            SPECIES &sp = Species[iptr->species];
            int nh = Species[iptr->species].nh;

            int *ivec = Atoms[gion].Qindex.data();

            int ncount = Atoms[gion].Qindex.size();

            for (int i=0; i < max_product * factor; i++)
                product[i] = 0.0;

            for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
            {

                /* Loop over states and accumulate charge */
                for (int istate = 0; istate < ct.num_states; istate++)
                {
                    double t1 = Kpts[kpt]->Kstates[istate].occupation[0] * Kpts[kpt]->kp.kweight;

                    for(int is = 0; is < ct.noncoll_factor; is++)
                    {    
                        for (int i = 0; i < ct.max_nl; i++)
                        {
                            sint[i + is * ct.max_nl] = Kpts[kpt]->newsint_local[(ct.noncoll_factor * istate + is)*num_nonloc_ions*ct.max_nl + ion * ct.max_nl + i];
                        }               /*end for i */
                    }

                    if(sp.is_spinorb) betapsi_so( (std::complex<double> *)sint, sint_tem, sp);

                    int idx = 0;
                    for (int i = 0; i < nh; i++)
                    {
                        for (int j = i; j < nh; j++)
                        {

                            double t2 = t1;
                            if(i != j) t2 = 2.0 * t1;

                            product[idx] += t2 * std::real(sint[i] * std::conj(sint[j]));

                            if(ct.noncoll)
                            {
                                std::complex<double> updown = 2.0 * sint[i] * std::conj(sint[j+ct.max_nl]);
                                product[idx + 1 * max_product] += t2 * std::real(updown);
                                product[idx + 2 * max_product] += t2 * std::imag(updown);
                                product[idx + 3 * max_product] += t2 * std::real(sint[i+ct.max_nl] * std::conj(sint[j+ct.max_nl]) );
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
                    for (int icount = 0; icount < ncount; icount++)
                    {
                        double Qr = Atoms[gion].augfunc[icount + idx * ncount];
                        augrho[ivec[icount]] += Qr * product[idx];
                        if(ct.noncoll)
                        {
                            augrho[ivec[icount]] += Qr * product[idx + 3 * max_product];
                            augrho[ivec[icount] + 1 * pbasis ] += Qr * product[idx + 1 * max_product];
                            augrho[ivec[icount] + 2 * pbasis ] += Qr * product[idx + 2 * max_product];
                            augrho[ivec[icount] + 3 * pbasis ] += Qr * (product[idx] - product[idx + 3 * max_product]);

                        }
                    }           /*end for icount */
                    idx++;
                }               /*end for j */
            }                   /*end for i */


        }                       /*end if */

    }                           /*end for ion */

    GlobalSums(augrho, pbasis * factor, pct.kpsub_comm);

    delete [] sint;
    delete [] product;

}

//  <beta|psi>  = sum_ fcoef_so * <beta|psi> 
// eq. 16 Corso and Conte, PRB 71, 115106 (2005)
void betapsi_so(std::complex<double> *sint, std::complex<double> *sint_tem, SPECIES &sp)
{
    for(int idx = 0; idx < 2*ct.max_nl; idx++) sint_tem[idx] = sint[idx];
    for(int idx = 0; idx < 2*ct.max_nl; idx++) sint[idx] = 0.0;

    for(int ih = 0; ih < sp.nh; ih++)
        for(int is = 0; is < 2; is++)
        {
            for(int jh = 0; jh < sp.nh; jh++)
                for(int is1 = 0; is1 < 2; is1++)
                {
                    sint[is * ct.max_nl + ih] += 
                        sp.fcoef_so[ih][jh][is*2+is1] * sint_tem[jh + is1 * ct.max_nl]; 
                }
        }

}

