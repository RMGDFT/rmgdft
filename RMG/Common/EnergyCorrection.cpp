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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "transition.h"
#include "Solvers.h"
#include "../Headers/prototypes.h"
#include "RmgParallelFft.h"
#include "RmgSumAll.h"


// Computes an approximate correction that helps make the total energy variational when using multigrid
// and potential acceleration to solve the Kohn-Sham equations. When using this combination as a solver
// the potential and the eigenvalues are updated continuously as the states are processed. This effectively
// means that there is a mixing of the input and output densities that is specfic to each state. In order
// to reduce this effect the initial values of the eigenvalues for each state are stored in MgEigState
// and used to generate a correction factor to be applied in GetTe.

template double EnergyCorrection<double> (Kpoint<double> **,
              double *, double *, double *, double *);
template double EnergyCorrection<std::complex<double> > (Kpoint<std::complex<double>> **,
              double *, double *, double *, double *);

template <typename OrbitalType> double EnergyCorrection (Kpoint<OrbitalType> **Kptr,
          double *rho, double *new_rho, double *vh, double *vh_in)

{

    int nspin = ct.spin_flag + 1;
    double ec = 0.0;
    bool potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));

    for (int is = 0; is < nspin; is++)
    {
        double ccharge = 0.0;
        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {

            Kpoint<OrbitalType> *kptr = Kptr[kpt];
            double t1 = 0.0;
            for (int st = 0; st < kptr->nstates; st++)
            {
                double mix = 1.0;
                if(potential_acceleration) mix = ccharge / ct.nel;
                t1 += (kptr->Kstates[st].occupation[is] * mix * (kptr->Kstates[st].feig[is] - kptr->Kstates[st].eig[is]));
                ccharge += kptr->Kstates[st].occupation[is];
            }
            ec += t1 * kptr->kweight;

        }
    }

    ec = RmgSumAll(ec, pct.kpsub_comm);
    return ec;
}


