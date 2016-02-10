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
#include "vdW.h"
#include "RmgException.h"
#include "xc.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "Pw.h"
#include "pfft.h"

#if USE_PFFT

// Used to performa a parallel interpolation from the wavefunction grid to the 
// potential grid using a phase shifting technique
void FftInterpolation (BaseGrid &G, double *coarse, double *fine, int ratio)
{

  int pbasis_c = G.get_P0_BASIS(1);

  if(ratio == 1) {
      for(int ix = 0;ix < pbasis_c;ix++) fine[ix] = coarse[ix];
      return;
  }

  ptrdiff_t n[3];
  n[0] = G.get_NX_GRID(1);
  n[1] = G.get_NY_GRID(1);
  n[2] = G.get_NZ_GRID(1);


  int dimx_c = G.get_PX0_GRID(1); 
  int dimy_c = G.get_PY0_GRID(1); 
  int dimz_c = G.get_PZ0_GRID(1); 
  int dimx_f = G.get_PX0_GRID(ratio); 
  int dimy_f = G.get_PY0_GRID(ratio); 
  int dimz_f = G.get_PZ0_GRID(ratio); 
  double scale = 1.0 / (double)(n[0]*n[1]*n[2]);
  double dratio = (double)ratio;

  // Array strides
  int incx_c = dimy_c * dimz_c;
  int incy_c = dimz_c;
  int incx_f = dimy_f * dimz_f;
  int incy_f = dimz_f;

  std::complex<double> *base_coarse = new std::complex<double>[pbasis_c];
  std::complex<double> *shifted_coarse = new std::complex<double>[pbasis_c];
  std::complex<double> *backshifted_coarse = new std::complex<double>[pbasis_c];
  pfft_plan forw, inv;

  forw = pfft_plan_dft_3d(n,
                          (double (*)[2])base_coarse,
                          (double (*)[2])base_coarse,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

  inv = pfft_plan_dft_3d(n,
                          (double (*)[2])shifted_coarse,
                          (double (*)[2])backshifted_coarse,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

  // Get offset of first grid point on this processor into global grid
  int offset_x, offset_y, offset_z;
  G.find_node_offsets(G.get_rank(), n[0], n[1], n[2], &offset_x, &offset_y, &offset_z);


  // Get the forward transform
  for(int ix = 0;ix < pbasis_c;ix++) base_coarse[ix] = std::complex<double>(coarse[ix], 0.0);
  pfft_execute_dft(forw, (double (*)[2])base_coarse, (double (*)[2])base_coarse);

  // Loop over phase shifts.
  for(int ix = 0;ix < ratio;ix++) {

      for(int iy = 0;iy < ratio;iy++) {

          for(int iz = 0;iz < ratio;iz++) {

              // Multiply coarse transform by phase shifts and 
              int idx = 0;
              for(int ixx = 0;ixx < dimx_c;ixx++) {

                  int p1 = offset_x + ixx;
                  if(p1 > n[0]/2) p1 -= n[0];
                  double rp1 = (double)(p1)*(double)(ix)/dratio / (double)n[0];

                  for(int iyy = 0;iyy < dimy_c;iyy++) {

                      int p2 = offset_y + iyy;
                      if(p2 > n[1]/2) p2 -= n[1];
                      double rp2 = (double)(p2)*(double)(iy)/dratio / (double)n[1];

                      for(int izz = 0;izz < dimz_c;izz++) {
                          int p3 = offset_z + izz;
                          if(p3 > n[2]/2) p3 -= n[2];
                          double rp3 = (double)p3*(double)(iz)/dratio / (double)n[2];
                          double theta = 2.0*PI*(rp1 + rp2 + rp3);
                          std::complex<double> phase = std::complex<double>(cos(theta), sin(theta));
                          shifted_coarse[idx] = phase * base_coarse[idx];
                          idx++;

                      }
                  }
              }

              // Backtransform phase shifted coefficients
              pfft_execute_dft(inv, (double (*)[2])shifted_coarse, (double (*)[2])backshifted_coarse);

              // Pack interpolated values into the fine grid
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      for(int izz = 0;izz < dimz_c;izz++) {
                          fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                             scale * (std::real(backshifted_coarse[ixx*incx_c + iyy*incy_c + izz]));

                      }
                  }
              }

          }

      }

  }



  pfft_destroy_plan(forw);
  pfft_destroy_plan(inv);
 
  delete [] backshifted_coarse;
  delete [] shifted_coarse;
  delete [] base_coarse;

}
#endif
