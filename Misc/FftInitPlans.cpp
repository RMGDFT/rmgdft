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

    int size = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    coarse_pwaves = new Pw(*Rmg_G, Rmg_L, 1, false);
    fine_pwaves = new Pw(*Rmg_G, Rmg_L, Rmg_G->default_FG_RATIO, false);

    std::complex<double> *tx = new std::complex<double>[size];
    grid[0] = Rmg_G->get_NX_GRID(1);
    grid[1] = Rmg_G->get_NY_GRID(1);
    grid[2] = Rmg_G->get_NZ_GRID(1);

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

}

#endif
