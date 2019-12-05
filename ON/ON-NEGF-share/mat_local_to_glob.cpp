/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include "LocalObject.h"
#include "blas.h"
#include "RmgException.h"
#include "transition.h"
#include "prototypes_on.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"


void mat_local_to_glob(double *mat_local, double *mat_glob, LocalObject<double> &A, LocalObject<double> &B, 
    int st_start1, int st_end1, int st_start2, int st_end2, bool reduce_flag)
{

//  st_start1 ... are global index for orbitals 
    int na = st_end1 - st_start1;
    int nb = st_end2 - st_start2;

    if(na < 1 || nb < 1) return;
    if(A.density != B.density)
        throw RmgFatalException() << "density is different "<< " at file " << __FILE__ << "\n";
        

    for(int idx = 0; idx < na * nb; idx++) mat_glob[idx] = 0.0;
    for(int st1 = st_start1; st1 < st_end1; st1++)
    {
        int st1_local = A.index_global_to_proj[st1];
        if (st1_local < 0 ) continue;

        for(int st2 = st_start2; st2 < st_end2; st2++)
        {
            int st2_local = B.index_global_to_proj[st2];
            if (st2_local < 0 ) continue;

            mat_glob[ (st2-st_start2) * na + st1-st_start1 ] += mat_local[st2_local * A.num_thispe + st1_local];
        }
    }

    if(reduce_flag)
    {
        int idx = na * nb;
        MPI_Allreduce(MPI_IN_PLACE, mat_glob, idx, MPI_DOUBLE, MPI_SUM, A.comm);
    }
}
