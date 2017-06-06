/*
 *
 * Copyright 2016 The RMG Project Developers. See the COPYRIGHT file 
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
#include "const.h"
#include "params.h"
#include "Atomic.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "blas.h"


#include "FiniteDiff.h"
#include "transition.h"


// Scf correction for forces as described by Chan, Bohnen and Ho
// PRB v47, #8

void CorrectForces (double * vh, double *vh_in, double *vxc, double *vxc_in, double *force)
{

    if(!ct.localize_localpp) pct.num_loc_ions = ct.num_ions;

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    double *gx = new double[3*FP0_BASIS];
    double *gy = gx + FP0_BASIS;
    double *gz = gy + FP0_BASIS;
    double *force_tmp = new double[pct.num_loc_ions * 3];
    double *dvh = new double[FP0_BASIS];

    int ithree = 3;
    double alpha = get_vel_f(), zero = 0.0;

    for(int ix = 0;ix < FP0_BASIS;ix++) dvh[ix] = vh_in[ix] + vxc_in[ix]  - vh[ix] - vxc[ix];

    ApplyGradient (dvh, gx, gy, gz, ct.force_grad_order, "Fine");

    double *dum_array = new double[FP0_BASIS];
    
    if(ct.localize_localpp)
        InitLocalObject (dum_array, pct.localatomicrho, ATOMIC_RHO, true);
    else
        InitDelocalizedObject (dum_array, pct.localatomicrho, ATOMIC_RHO, true);

    dgemm("T", "N", &ithree, &pct.num_loc_ions, &FP0_BASIS, &alpha, gx, &FP0_BASIS, 
            pct.localatomicrho, &FP0_BASIS, &zero, force_tmp, &ithree); 

    delete [] pct.localatomicrho;
    delete [] dum_array;

    for(int ion1 = 0; ion1 <pct.num_loc_ions; ion1++)
    {
        int ion = pct.loc_ions_list[ion1];

        force[ion *3 + 0] = force_tmp[ion1 *3 + 0];
        force[ion *3 + 1] = force_tmp[ion1 *3 + 1];
        force[ion *3 + 2] = force_tmp[ion1 *3 + 2];

    }


    delete [] dvh;
    delete [] force_tmp;
    delete [] gx;


}

