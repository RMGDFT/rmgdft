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


#include <stdio.h>
#include <stdlib.h>
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
#include "Kpoint.h"
#include "GpuAlloc.h"
#include "RmgGemm.h"
#include "blas.h"



template void Betaxpsi<double>(Kpoint<double> *, int, int, double *);
template void Betaxpsi<std::complex<double> >(Kpoint<std::complex<double>> *, int, int, std::complex<double> *);


template <typename KpointType>
void Betaxpsi (Kpoint<KpointType> *kptr, int first_state, int nstates, KpointType *sint_local)
{
    kptr->BetaProjector->project(kptr, sint_local, first_state, nstates, kptr->nl_weight, kptr->pbasis);
}
