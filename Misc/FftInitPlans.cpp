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
#include "RmgParallelFft.h"

// Declared here with extern declarations in transition.h
Pw *coarse_pwaves, *fine_pwaves, *ewald_pwaves;

// Initializes common plans and plane wave objects for reuse.
void FftInitPlans(void)
{
    coarse_pwaves = new Pw(*Rmg_G, Rmg_L, 1, false);
    fine_pwaves = new Pw(*Rmg_G, Rmg_L, Rmg_G->default_FG_RATIO, false);

    if(Rmg_G->default_FG_RATIO > 1)
    {
        ewald_pwaves = fine_pwaves;
    }
    else
    {
        ewald_pwaves = new Pw(*Rmg_G, Rmg_L, 2, false);
    }

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

