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

void LO_x_LO(LocalObject<double> &A, LocalObject<double> &B, double *mat, BaseGrid &Rmg_G)
{
    double *mat_local;

    int na = A.num_thispe ;
    int nb = B.num_thispe;
    if(A.density != B.density)
        throw RmgFatalException() << "density is different "<< " at file " << __FILE__ << "\n";
        
    int P0_BASIS = Rmg_G.get_P0_BASIS(A.density);
    
    mat_local = new double[na*nb];
    double one = 1.0, zero = 0.0;

    dgemm("T", "N", &na, &nb, &P0_BASIS, &one, A.storage_proj, &P0_BASIS,
                B.storage_proj, &P0_BASIS, &zero, mat_local, &na);

    for (int j = 0; j < B.num_tot; j ++)
    for (int i = 0; i < A.num_tot; i ++)
        mat[j * A.num_tot + i] = 0.0;
    
    double t1 = Rmg_G.get_NX_GRID(A.density); 
    t1 *= Rmg_G.get_NY_GRID(A.density); 
    t1 *= Rmg_G.get_NZ_GRID(A.density); 

    double vol = Rmg_L.get_omega() /t1;
    for (int j = 0; j < B.num_thispe; j++)
    for (int i = 0; i < A.num_thispe; i++)
    {
        int i_glob = A.index_proj_to_global[i];
        int j_glob = B.index_proj_to_global[j];
        mat[j_glob * A.num_tot + i_glob] = mat_local[j*A.num_thispe+i] * vol;

    }

    int idx = A.num_tot * B.num_tot;
    MPI_Allreduce(MPI_IN_PLACE, (double *)mat, idx, MPI_DOUBLE, MPI_SUM, A.comm);
    delete [] mat_local;

}

