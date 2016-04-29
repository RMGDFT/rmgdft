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
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "transition.h"

#if USE_PFFT
#include "RmgParallelFft.h"

// Declared here with extern declarations in transition.h
pfft_plan forward_coarse, backward_coarse, forward_fine, backward_fine;
Pw *coarse_pwaves, *fine_pwaves;


// Initializes common plans and plane wave objects for reuse.
void FftInitPlans(void)
{

    ptrdiff_t grid[3]; 
    ptrdiff_t local_ni[3], local_i_start[3], local_no[3], local_o_start[3];
    int np[3];

    // First we check if standard decompositions will work
    np[0] = Rmg_G->get_PE_X();
    np[1] = Rmg_G->get_PE_Y();
    np[2] = Rmg_G->get_PE_Z();
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);

    pfft_init();
    if( pfft_create_procmesh(3, pct.grid_comm, np, &pct.pfft_comm) ) {
        RmgFatalException() << "Problem initializing PFFT in " << __FILE__ << " at line " << __LINE__ << ".\n";
    }

    grid[0] = Rmg_G->get_NX_GRID(1);
    grid[1] = Rmg_G->get_NY_GRID(1);
    grid[2] = Rmg_G->get_NZ_GRID(1);

    // get array sizes
    int local_size = pfft_local_size_dft_3d(grid, pct.pfft_comm, PFFT_TRANSPOSED_NONE,
                                          local_ni, local_i_start, local_no, local_o_start);

//printf("local_size = %d\n", local_size);
//printf("local_ni = %d %d  %d  %d  %d  %d  %d\n", local_size, local_ni[0], local_ni[1], local_ni[2], dimx, dimy, dimz);
//printf("local_ni_start = %d  %d  %d\n", local_i_start[0], local_i_start[1], local_i_start[2]);
//printf("local_no = %d  %d  %d  %d\n", local_size, local_no[0], local_no[1], local_no[2]);
//printf("local_no_start = %d  %d  %d\n", local_o_start[0], local_o_start[1], local_o_start[2]);


    std::vector<int> zfactors = {1};
    GetPrimeFactors(zfactors, np[2], np[2]);


    coarse_pwaves = new Pw(*Rmg_G, Rmg_L, 1, false, pct.pfft_comm);
    fine_pwaves = new Pw(*Rmg_G, Rmg_L, Rmg_G->default_FG_RATIO, false, pct.pfft_comm);

    // Fine grid size is large enough for both coarse and fine
    int size = fine_pwaves->local_size;
//printf("DDDDDD %d\n",size);
    std::complex<double> *tx = new std::complex<double>[size];

    forward_coarse =  pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    backward_coarse = pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    grid[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    grid[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    grid[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);

    forward_fine =  pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    backward_fine = pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    delete [] tx;
}

#endif
