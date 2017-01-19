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
#include "remap_3d.h"
#include "RmgParallelFft.h"

// Declared here with extern declarations in transition.h
pfft_plan forward_coarse, backward_coarse, forward_fine, backward_fine;
Pw *coarse_pwaves, *fine_pwaves;

// Initializes common plans and plane wave objects for reuse.
void FftInitPlans(void)
{

    ptrdiff_t grid[3]; 
    ptrdiff_t coarse_ni[3], coarse_i_start[3], coarse_no[3], coarse_o_start[3];
    ptrdiff_t fine_ni[3], fine_i_start[3], fine_no[3], fine_o_start[3];
    int np[3];

    // First we check if standard decomposition with pfft will work
    np[0] = Rmg_G->get_PE_X();
    np[1] = Rmg_G->get_PE_Y();
    np[2] = Rmg_G->get_PE_Z();


    pfft_init();
    if( pfft_create_procmesh(3, pct.grid_comm, np, &pct.pfft_comm) ) {
        RmgFatalException() << "Problem initializing PFFT in " << __FILE__ << " at line " << __LINE__ << ".\n";
    }

    // See if we can use pfft without remapping. In and out arrays must be equal in size.
    grid[0] = Rmg_G->get_NX_GRID(1);
    grid[1] = Rmg_G->get_NY_GRID(1);
    grid[2] = Rmg_G->get_NZ_GRID(1);
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);

    int coarse_size = pfft_local_size_dft_3d(grid, pct.pfft_comm, PFFT_TRANSPOSED_NONE,
                                          coarse_ni, coarse_i_start, coarse_no, coarse_o_start);

    int coarse_remap = (coarse_ni[0] != dimx) || (coarse_ni[1] != dimy) || (coarse_ni[2] != dimz) ||
                (coarse_no[0] != dimx) || (coarse_no[1] != dimy) || (coarse_no[2] != dimz);
    MPI_Allreduce(MPI_IN_PLACE, &coarse_remap, 1, MPI_INT, MPI_SUM, pct.grid_comm);

    coarse_pwaves = new Pw(*Rmg_G, Rmg_L, 1, false, pct.pfft_comm);
    coarse_pwaves->remap_local_size = coarse_size;
    coarse_pwaves->fwd_remap = NULL;
    coarse_pwaves->inv_remap = NULL;

    // Create remap plans if needed
    int pxoffset, pyoffset, pzoffset;
    if(coarse_remap) {
        if(pct.gridpe == 0) printf("Remapping coarse grids for parallel fft.\n");
        
        Rmg_G->find_node_offsets(pct.gridpe, grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);    
        // We only treat the double precision complex case
        coarse_pwaves->fwd_remap = remap_3d_create_plan(
                           pct.grid_comm,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           coarse_i_start[2], coarse_i_start[2] + coarse_ni[2] - 1,
                           coarse_i_start[1], coarse_i_start[1] + coarse_ni[1] - 1,
                           coarse_i_start[0], coarse_i_start[0] + coarse_ni[0] - 1,
                           sizeof(std::complex<double>)/sizeof(double), 0, 1, 2);

        coarse_pwaves->inv_remap = remap_3d_create_plan(
                           pct.grid_comm,
                           coarse_o_start[2], coarse_o_start[2] + coarse_no[2] - 1,
                           coarse_o_start[1], coarse_o_start[1] + coarse_no[1] - 1,
                           coarse_o_start[0], coarse_o_start[0] + coarse_no[0] - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           sizeof(std::complex<double>)/sizeof(double), 0, 1, 2);

    }



    grid[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    grid[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    grid[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
    dimx = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    dimy = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    dimz = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int fine_size = pfft_local_size_dft_3d(grid, pct.pfft_comm, PFFT_TRANSPOSED_NONE,
                                          fine_ni, fine_i_start, fine_no, fine_o_start);

    int fine_remap = (fine_ni[0] != dimx) || (fine_ni[1] != dimy) || (fine_ni[2] != dimz) ||
                (fine_no[0] != dimx) || (fine_no[1] != dimy) || (fine_no[2] != dimz);
    MPI_Allreduce(MPI_IN_PLACE, &fine_remap, 1, MPI_INT, MPI_SUM, pct.grid_comm);

    fine_pwaves = new Pw(*Rmg_G, Rmg_L, Rmg_G->default_FG_RATIO, false, pct.pfft_comm);
    fine_pwaves->remap_local_size = fine_size;
    fine_pwaves->fwd_remap = NULL;
    fine_pwaves->inv_remap = NULL;

    if(fine_remap) {
        
        if(pct.gridpe == 0) printf("Remapping fine grids for parallel fft.\n");
        Rmg_G->find_node_offsets(pct.gridpe, grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);    
        // We only treat the double precision complex case
        fine_pwaves->fwd_remap = remap_3d_create_plan(
                           pct.grid_comm,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           fine_i_start[2], fine_i_start[2] + fine_ni[2] - 1,
                           fine_i_start[1], fine_i_start[1] + fine_ni[1] - 1,
                           fine_i_start[0], fine_i_start[0] + fine_ni[0] - 1,
                           sizeof(std::complex<double>)/sizeof(double), 0, 1, 2);

        fine_pwaves->inv_remap = remap_3d_create_plan(
                           pct.grid_comm,
                           fine_o_start[2], fine_o_start[2] + fine_no[2] - 1,
                           fine_o_start[1], fine_o_start[1] + fine_no[1] - 1,
                           fine_o_start[0], fine_o_start[0] + fine_no[0] - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           sizeof(std::complex<double>)/sizeof(double), 0, 1, 2);

    }


    // Fine grid size is large enough for both coarse and fine
    std::complex<double> *tx = new std::complex<double>[fine_size];

    grid[0] = Rmg_G->get_NX_GRID(1);
    grid[1] = Rmg_G->get_NY_GRID(1);
    grid[2] = Rmg_G->get_NZ_GRID(1);
    forward_coarse =  pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
    coarse_pwaves->forward_plan = &forward_coarse;

    backward_coarse = pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
    coarse_pwaves->backward_plan = &backward_coarse;


    grid[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    grid[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    grid[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
    forward_fine =  pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
    fine_pwaves->forward_plan = &forward_fine;

    backward_fine = pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
    fine_pwaves->backward_plan = &backward_fine;

    delete [] tx;

}

