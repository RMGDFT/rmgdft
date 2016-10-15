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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#if !(defined(_WIN32) || defined(_WIN64))
#include <unistd.h>
#else
#include <io.h>
#endif
#include <time.h>
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"
#include "RmgParallelFft.h"


template void WriteBGW ( double *, double *, double *, Kpoint<double> **);
template void WriteBGW ( double *, double *, double *, Kpoint<std::complex<double> > **);

    template <typename KpointType>
void WriteBGW ( double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr)
{
    WriteRhog(rho, rho_oppo);
    WriteVxcEig(vxc, Kptr);

}

