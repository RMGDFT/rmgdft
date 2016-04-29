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
#include <math.h>
#include <float.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "transition.h"

#if USE_PFFT
#include "RmgParallelFft.h"
#endif

template void PfftForward<double, std::complex<double> >(double *, std::complex<double> *, char *);

// This function performs a parallel forward FFT using the pfft library. If the actual
// processor grid is not directly compatible with the pfft requirements it handles the
// required data movement in order to perform the FFT.
template <typename InType, typename OutType>
void PfftForward (InType * in, OutType * out, char *density)
{
}
