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
#include "RmgSumAll.h"
#include "transition.h"
#include "GlobalSums.h"
#include "RmgParallelFft.h"

// On input performs a dft of x which is an array distributed in real space across all nodes
// using the plane wave structure defined in pwaves. It then computes a frequency histogram
// of all of the coefficients of the transform. This is intended as a debugging function
// and should not be called from a threaded region.
void FftFreqBin(double *x,   // IN:OUT  Input array in real space. Distributed across all nodes.
               Pw &pwaves,  // IN:     Plane wave structure that corresponds to the reciprocal space grid for x
               double *bins)  // IN:OUT   Allocated by calling array. Stores frequency bins. Calling array is
                              // respsponsible for making sure the array is big enough ((int)rint(pwaves.gmax) + 1)
 
{
#if 0
  ptrdiff_t grid[3];
  pfft_plan forw;
  grid[0] = pwaves.global_dimx;
  grid[1] = pwaves.global_dimy;
  grid[2] = pwaves.global_dimz;

  int pbasis = pwaves.pbasis;
  int nvecs = (int)rint(sqrt(pwaves.gmax)) + 1;
  int *binsize = new int[nvecs+1]();

  for(int i=0;i < nvecs;i++)bins[i] = 0.0;

  std::complex<double> *cvec = new std::complex<double>[pbasis]();
  if(&pwaves == coarse_pwaves) {
      forw = forward_coarse;
  }
  else if(&pwaves == fine_pwaves) {
      forw = forward_fine;
  }
  else {
      forw = pfft_plan_dft_3d(grid,
                                  (double (*)[2])cvec,
                                  (double (*)[2])cvec,
                                  pct.pfft_comm,
                                  PFFT_FORWARD,
                                  PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
  }

  for(int i = 0;i < pbasis;i++) cvec[i] = std::complex<double>(x[i], 0.0);
  pfft_execute_dft(forw, (double (*)[2])cvec, (double (*)[2])cvec);

  for(int ig=0;ig < pbasis;ig++) {
      double t1 = sqrt(std::norm(cvec[ig]));
      int gidx = (int)rint(sqrt(pwaves.gmags[ig]));
      if(gidx <= nvecs) {
          bins[gidx] += t1;
          binsize[gidx]++;
      }
  }

  MPI_Allreduce(MPI_IN_PLACE, bins, nvecs, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
  MPI_Allreduce(MPI_IN_PLACE, binsize, nvecs, MPI_INT, MPI_SUM, pct.grid_comm);

  for(int i = 0;i < nvecs;i++) {
      if(binsize[i]) bins[i] /= (double)binsize[i];
  }

  if((&pwaves != coarse_pwaves) && (&pwaves != fine_pwaves)) {
      pfft_destroy_plan(forw);
  }

  delete [] binsize;
  delete [] cvec;
#endif
}

