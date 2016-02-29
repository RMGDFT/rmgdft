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
#include "RmgSumAll.h"
#include "transition.h"

#if USE_PFFT

#include "RmgParallelFft.h"

// On input performs a dft of x which is an array distributed in real space across all node
// using the plane wave structure defined in pwaves. This is then filtered in g-space with
// the coefficients of all plane waves with a cutoff greater than the maximum value on the grid
// set to zero.
void FftFilter(double *x,   // IN:OUT  Input array in real space. Distributed across all nodes.
               Pw &pwaves,  // IN:     Plane wave structure that corresponds to the reciprocal space grid for x
               double factor)  // IN:  Plane wave filtering factor between (0.0, 1.). Closely corresponds to 
                            // maximum frequence on a grid of equivalent density.
{

  if((factor <= 0.0) || (factor > 1.0))
      throw RmgFatalException() << "FFT filtering factor must be between 0.0 and 1.0 " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

  ptrdiff_t densgrid[3];
  pfft_plan plan_forward, plan_back;
  densgrid[0] = pwaves.global_dimx;
  densgrid[1] = pwaves.global_dimy;
  densgrid[2] = pwaves.global_dimz;

  double g2cut = factor*factor*(sqrt(pwaves.gmax))*(sqrt(pwaves.gmax));
  int global_basis = pwaves.global_basis;
  int pbasis = pwaves.pbasis;

  std::complex<double> *crho = new std::complex<double>[pbasis];
  plan_forward = pfft_plan_dft_3d(densgrid,
                                  (double (*)[2])crho,
                                  (double (*)[2])crho,
                                  pct.pfft_comm,
                                  PFFT_FORWARD,
                                  PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

  plan_back = pfft_plan_dft_3d(densgrid,
                               (pfft_complex *)crho,
                               (pfft_complex *)crho,
                               pct.pfft_comm,
                               PFFT_BACKWARD,
                               PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

  for(int i = 0;i < pbasis;i++) crho[i] = std::complex<double>(x[i], 0.0);
  pfft_execute_dft(plan_forward, (double (*)[2])crho, (double (*)[2])crho);
  for(int ig=0;ig < pbasis;ig++) {
      if(pwaves.gmags[ig] > g2cut) {
          crho[ig] = std::complex<double>(0.0, 0.0);
      }
  }
  pfft_execute_dft(plan_back, (double (*)[2])crho, (double (*)[2])crho);
  for(int i = 0;i < pbasis;i++) x[i] = std::real(crho[i])/(double)global_basis;

  pfft_destroy_plan(plan_back);
  pfft_destroy_plan(plan_forward);
  delete [] crho;

}

#endif
