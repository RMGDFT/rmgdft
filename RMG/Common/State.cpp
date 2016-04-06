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
#include <complex>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "RmgSumAll.h"
#include <complex>
#include "../Headers/prototypes.h"


template State<double>::State(void);
template State<std::complex<double> >::State(void);

template void State<double>::normalize(double *, int );
template void State<std::complex <double> >::normalize(std::complex <double> *, int );
template void State<double>::set_storage(double *storage);
template void State<std::complex <double> >::set_storage(std::complex <double> *tpsi);
template bool State<double>::is_occupied(void);
template bool State<std::complex <double> >::is_occupied(void);

template <class StateType> State<StateType>::State(void)
{
    this->occupation[0] = 0.0;
    this->occupation[1] = 0.0;

}

template <class StateType> void State<StateType>::set_storage(StateType *storage)
{
    this->psi = storage;
}

template <class StateType> bool State<StateType>::is_occupied(void)
{
    if(!ct.spin_flag) return (this->occupation[0] > 0.0);
    return ((this->occupation[0] > 0.0) || (this->occupation[1] > 0.0));
}

template <class StateType> void State<StateType>::normalize(StateType *tpsi, int istate)
{
    double vel = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    vel = Rmg_L.get_omega() / vel;

    if(ct.norm_conserving_pp) {

        double sum1 = 0.0, sum2;
        for(int idx = 0;idx < this->Kptr->pbasis;idx++) {
            sum1 = sum1 + std::norm(tpsi[idx]);
        }

        sum2 = vel * RmgSumAll<double>(sum1, this->Kptr->comm);
        sum2 = sqrt(1.0 / sum2);

        StateType tscale(sum2);
        for(int idx = 0;idx < this->Kptr->pbasis;idx++) {
            tpsi[idx] = tpsi[idx] * tscale;
        }

    }
    else {

        int ion, nh, i, j, inh, sidx_local, nidx, oion;
        double sumbeta, sumpsi, sum, t1;
        double *qqq;
        StateType *sint;
        ION *iptr;

        sidx_local = istate * pct.num_nonloc_ions * ct.max_nl;

        sumpsi = 0.0;
        sumbeta = 0.0;

        nidx = -1;
        for (ion = 0; ion < pct.num_owned_ions; ion++)
        {

            oion = pct.owned_ions_list[ion];
            
            iptr = &ct.ions[oion];
           
            nh = ct.sp[iptr->species].nh;
            
            /* Figure out index of owned ion in nonloc_ions_list array*/
            do {
                
                nidx++;
                if (nidx >= pct.num_nonloc_ions)
                    //error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
                    rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in pct.nonloc_ions_list");
            
            } while (pct.nonloc_ions_list[nidx] != oion);

            qqq = pct.qqq[oion];


            /*nidx adds offset due to current ion*/
            sint = &this->Kptr->newsint_local[sidx_local + nidx * ct.max_nl];

            for (i = 0; i < nh; i++)
            {
                inh = i * nh;
                for (j = 0; j < nh; j++)
                {
                    if (qqq[inh + j] != 0.0)
                    {

                        sumbeta += qqq[inh + j] * std::real((std::real(sint[i]) * std::real(sint[j]) + std::imag(sint[i]) * std::imag(sint[j])));
                        sumbeta += qqq[inh + j] * std::real((std::real(sint[i]) * std::imag(sint[j]) - std::imag(sint[i]) * std::real(sint[j])));

                    }
                }
            }
        }


        for (int idx = 0; idx < this->Kptr->pbasis; idx++)
        {
            sumpsi += std::norm(tpsi[idx]);
        }

        sum = RmgSumAll<double>(vel * sumpsi + sumbeta, this->Kptr->comm);
        sum = 1.0 / sum;

        if (sum < 0.0)
        {
            rmg_printf ("the %dth state is wrong\n", istate);
            rmg_error_handler (__FILE__, __LINE__, "<psi|S|psi> cann't be negative");
        }

        t1 = sqrt (sum);
        for(int idx = 0;idx < this->Kptr->pbasis;idx++) {
            tpsi[idx] = tpsi[idx] * t1;
        }

        /* update <beta|psi> - Local version*/
        
        for (ion = 0; ion < pct.num_nonloc_ions; ion++)
        {

            StateType *tsint_ptr = &this->Kptr->newsint_local[sidx_local + ion * ct.max_nl];
            for(int jdx = 0;jdx < ct.max_nl;jdx++) {
                tsint_ptr[jdx] = tsint_ptr[jdx] * t1;
            }

        }

    }
}

