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
#include "transition.h"
#include "GlobalSums.h"
#include "RmgException.h"
#include "fft3d.h"
#include "RmgParallelFft.h"

// Declared here with extern declarations in transition.h
Pw *coarse_pwaves, *fine_pwaves;
struct fft_plan_3d *fft_forward_coarse, *fft_backward_coarse, *fft_forward_fine, *fft_backward_fine;

// Initializes common plans and plane wave objects for reuse.
void FftInitPlans(void)
{

    int grid[3]; 
    pct.pfft_comm = pct.grid_comm;

    // See if we can use pfft without remapping. In and out arrays must be equal in size.
    grid[0] = Rmg_G->get_NX_GRID(1);
    grid[1] = Rmg_G->get_NY_GRID(1);
    grid[2] = Rmg_G->get_NZ_GRID(1);
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);
   
    coarse_pwaves = new Pw(*Rmg_G, Rmg_L, 1, false, pct.grid_comm);

    int pxoffset, pyoffset, pzoffset, nbuf, scaled=false, permute=0, usecollective=true;
    Rmg_G->find_node_offsets(pct.gridpe, grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);    
    fft_forward_coarse = fft_3d_create_plan(pct.grid_comm,
                           grid[2], grid[1], grid[0],
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           scaled, permute, &nbuf, usecollective);

    fft_backward_coarse = fft_forward_coarse;
    coarse_pwaves->fft_forward_plan = fft_forward_coarse;
    coarse_pwaves->fft_backward_plan = fft_forward_coarse;


    grid[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    grid[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    grid[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
    dimx = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    dimy = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    dimz = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);

    fine_pwaves = new Pw(*Rmg_G, Rmg_L, Rmg_G->default_FG_RATIO, false, pct.grid_comm);

    Rmg_G->find_node_offsets(pct.gridpe, grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);    
    fft_forward_fine = fft_3d_create_plan(pct.grid_comm,
                           grid[2], grid[1], grid[0],
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           scaled, permute, &nbuf, usecollective);

    fft_backward_fine = fft_forward_fine;
    fine_pwaves->fft_forward_plan = fft_forward_fine;
    fine_pwaves->fft_backward_plan = fft_forward_fine;

}

/* ----------------------------------------------------------------------
   Create plan for performing a 3d FFT

   Arguments:
   comm                 MPI communicator for the P procs which own the data
   nfast,nmid,nslow     size of global 3d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute once = mid->fast, slow->mid, fast->slow
                          2 = permute twice = slow->fast, fast->mid, mid->slow
   nbuf                 returns size of internal storage buffers used by FFT
   usecollective        use collective MPI operations for remapping data
------------------------------------------------------------------------- */

