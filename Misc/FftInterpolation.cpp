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
#include "RmgParallelFft.h"
#include <boost/math/special_functions/powm1.hpp>

void Fftpack_coarse_to_fine(std::complex<double> *coarse, double *fine,
                            int dimx_c, int dimy_c, int dimz_c, 
                            int xshift, int yshift, int zshift, 
                            double scale, int ratio);

// Used to performa a parallel interpolation from the wavefunction grid to the 
// potential grid using a phase shifting technique
void FftInterpolation (BaseGrid &G, double *coarse, double *fine, int ratio, bool use_sqrt)
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
  int dimy_f = G.get_PY0_GRID(ratio); 
  int dimz_f = G.get_PZ0_GRID(ratio); 
  double scale = 1.0 / (double)(n[0]*n[1]*n[2]);
  double dratio = (double)ratio;
  double rootrho;

  // Array strides
  int incx_c = dimy_c * dimz_c;
  int incy_c = dimz_c;
  int incx_f = dimy_f * dimz_f;
  int incy_f = dimz_f;

  std::complex<double> *base_coarse = new std::complex<double>[pbasis_c];
  std::complex<double> *shifted_coarse = new std::complex<double>[pbasis_c];
  std::complex<double> *backshifted_coarse = new std::complex<double>[pbasis_c];

  // Get offset of first grid point on this processor into global grid
  int offset_x, offset_y, offset_z;
  G.find_node_offsets(G.get_rank(), n[0], n[1], n[2], &offset_x, &offset_y, &offset_z);

  // Get the forward transform
  if(use_sqrt)
  {
      for(int ix = 0;ix < pbasis_c;ix++) base_coarse[ix] = std::complex<double>(sqrt(coarse[ix]), 0.0);
  }
  else
  {
      for(int ix = 0;ix < pbasis_c;ix++) base_coarse[ix] = std::complex<double>(coarse[ix], 0.0);
  }

  coarse_pwaves->FftForward(base_coarse, base_coarse);

  // Loop over phase shifts.
  for(int ix = 0;ix < ratio;ix++) {

      for(int iy = 0;iy < ratio;iy++) {

          for(int iz = 0;iz < ratio;iz++) {

              // Multiply coarse transform by phase shifts and 
              int idx = 0;
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  int p1 = offset_x + ixx;
                  if(p1 > n[0]/2) p1 -= n[0];
                  //if(ct.sqrt_interpolation && (p1 == n[0]/2)) p1 = -p1;
                  double rp1 = (double)(p1)*(double)(ix)/dratio / (double)n[0];
                  double theta_x = 2.0*PI*rp1;
                  std::complex<double> phase_x = std::complex<double>(cos(theta_x), sin(theta_x));
//                  if(p1 == n[0]/2) phase_x = std::real(phase_x);
                  //if(p1 == n[0]/2) phase_x = 1.0;

                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      int p2 = offset_y + iyy;
                      if(p2 > n[1]/2) p2 -= n[1];
                      //if(ct.sqrt_interpolation && (p2 == n[1]/2)) p2 = -p2;
                      double rp2 = (double)(p2)*(double)(iy)/dratio / (double)n[1];
                      double theta_y = 2.0*PI*rp2;
                      std::complex<double> phase_y = std::complex<double>(cos(theta_y), sin(theta_y));
                      //if(p2 == n[1]/2) phase_y = std::real(phase_y);
                      //if(p2 == n[1]/2) phase_y = 1.0;

                      for(int izz = 0;izz < dimz_c;izz++) {
                          int p3 = offset_z + izz;
                          if(p3 > n[2]/2) p3 -= n[2];
                          //if(ct.sqrt_interpolation && (p3 == n[2]/2)) p3 = -p3;
                          double rp3 = (double)p3*(double)(iz)/dratio / (double)n[2];
                          double theta_z = 2.0*PI*rp3;
                          std::complex<double> phase_z = std::complex<double>(cos(theta_z), sin(theta_z));
                          //if(p3 == n[2]/2) phase_z = std::real(phase_z);
                          //if(p3 == n[2]/2) phase_z = 1.0;
                          std::complex<double> phase = phase_x * phase_y * phase_z;
                          shifted_coarse[idx] = phase * base_coarse[idx];
                          idx++;

                      }
                  }
              }

              // Backtransform phase shifted coefficients
              coarse_pwaves->FftInverse(shifted_coarse, backshifted_coarse);

              // Pack interpolated values into the fine grid
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      for(int izz = 0;izz < dimz_c;izz++) {
                          if(use_sqrt)
                          {
                              rootrho =  std::abs(backshifted_coarse[ixx*incx_c + iyy*incy_c + izz]);
                              fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                                  (scale*scale)*(rootrho*rootrho);
                          }
                          else
                          {
                              rootrho =  std::real(backshifted_coarse[ixx*incx_c + iyy*incy_c + izz]);
                              fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                                  scale*rootrho;
                          }

                      }
                  }
              }

          }

      }

  }

  // Loop over phase shifts.
  for(int ix = 0;ix < ratio;ix++) {

      for(int iy = 0;iy < ratio;iy++) {

          for(int iz = 0;iz < ratio;iz++) {

              // Multiply coarse transform by phase shifts and 
              int idx = 0;
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  int p1 = offset_x + ixx;
                  if(p1 > n[0]/2) p1 -= n[0];
                  if(ct.sqrt_interpolation && (p1 == n[0]/2)) p1 = -p1;
                  double rp1 = (double)(p1)*(double)(ix)/dratio / (double)n[0];
                  double theta_x = 2.0*PI*rp1;
                  std::complex<double> phase_x = std::complex<double>(cos(theta_x), sin(theta_x));
//                  if(p1 == n[0]/2) phase_x = std::real(phase_x);
                  //if(p1 == n[0]/2) phase_x = 1.0;

                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      int p2 = offset_y + iyy;
                      if(p2 > n[1]/2) p2 -= n[1];
                     if(ct.sqrt_interpolation && (p2 == n[1]/2)) p2 = -p2;
                      double rp2 = (double)(p2)*(double)(iy)/dratio / (double)n[1];
                      double theta_y = 2.0*PI*rp2;
                      std::complex<double> phase_y = std::complex<double>(cos(theta_y), sin(theta_y));
                      //if(p2 == n[1]/2) phase_y = std::real(phase_y);
                      //if(p2 == n[1]/2) phase_y = 1.0;

                      for(int izz = 0;izz < dimz_c;izz++) {
                          int p3 = offset_z + izz;
                          if(p3 > n[2]/2) p3 -= n[2];
                         if(ct.sqrt_interpolation && (p3 == n[2]/2)) p3 = -p3;
                          double rp3 = (double)p3*(double)(iz)/dratio / (double)n[2];
                          double theta_z = 2.0*PI*rp3;
                          std::complex<double> phase_z = std::complex<double>(cos(theta_z), sin(theta_z));
                          //if(p3 == n[2]/2) phase_z = std::real(phase_z);
                          //if(p3 == n[2]/2) phase_z = 1.0;
                          std::complex<double> phase = phase_x * phase_y * phase_z;
                          shifted_coarse[idx] = phase * base_coarse[idx];
                          idx++;

                      }
                  }
              }

              // Backtransform phase shifted coefficients
              coarse_pwaves->FftInverse(shifted_coarse, backshifted_coarse);

              // Pack interpolated values into the fine grid
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      for(int izz = 0;izz < dimz_c;izz++) {
                          if(use_sqrt)
                          {
                              rootrho =  std::abs(backshifted_coarse[ixx*incx_c + iyy*incy_c + izz]);
                              fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                              (fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] + 
                                  (scale*scale)*(rootrho*rootrho))/2.0;
                          }
                          else
                          {
                              rootrho =  std::real(backshifted_coarse[ixx*incx_c + iyy*incy_c + izz]);
                              fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                                  scale*rootrho;
                          }

                      }
                  }
              }

          }

      }

  }


  delete [] backshifted_coarse;
  delete [] shifted_coarse;
  delete [] base_coarse;

}

void Fftpack_coarse_to_fine(std::complex<double> *coarse, double *fine, 
        int dimx_c, int dimy_c, int dimz_c, 
        int xshift, int yshift, int zshift, 
        double scale, int ratio)
{

    // Array strides
    int incx_c = dimy_c * dimz_c;
    int incy_c = dimz_c;

    int incx_f = ratio * ratio * dimy_c * dimz_c;
    int incy_f = ratio * dimz_c;

    // Pack interpolated values into the fine grid
    for(int ixx = 0;ixx < dimx_c;ixx++) {
        for(int iyy = 0;iyy < dimy_c;iyy++) {
            for(int izz = 0;izz < dimz_c;izz++) {
                fine[(ratio*ixx + xshift)*incx_f + (ratio*iyy + yshift)*incy_f + (ratio*izz + zshift)] = 
                    scale * (std::real(coarse[ixx*incx_c + iyy*incy_c + izz]));

            }
        }
    }
}

void FftInterpolation (BaseGrid &G, std::complex<double> *coarse, std::complex<double> *fine, int ratio, bool use_sqrt)
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

  // Get offset of first grid point on this processor into global grid
  int offset_x, offset_y, offset_z;
  G.find_node_offsets(G.get_rank(), n[0], n[1], n[2], &offset_x, &offset_y, &offset_z);

  // Get the forward transform
  for(int ix = 0;ix < pbasis_c;ix++) base_coarse[ix] = coarse[ix];

  coarse_pwaves->FftForward(base_coarse, base_coarse);

  // Loop over phase shifts.
  for(int ix = 0;ix < ratio;ix++) {

      for(int iy = 0;iy < ratio;iy++) {

          for(int iz = 0;iz < ratio;iz++) {

              // Multiply coarse transform by phase shifts and 
              int idx = 0;
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  int p1 = offset_x + ixx;
                  if(p1 > n[0]/2) p1 -= n[0];
                  //if(ct.sqrt_interpolation && (p1 == n[0]/2)) p1 = -p1;
                  double rp1 = (double)(p1)*(double)(ix)/dratio / (double)n[0];
                  double theta_x = 2.0*PI*rp1;
                  std::complex<double> phase_x = std::complex<double>(cos(theta_x), sin(theta_x));
//                  if(p1 == n[0]/2) phase_x = std::real(phase_x);
                  //if(p1 == n[0]/2) phase_x = 1.0;

                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      int p2 = offset_y + iyy;
                      if(p2 > n[1]/2) p2 -= n[1];
                      //if(ct.sqrt_interpolation && (p2 == n[1]/2)) p2 = -p2;
                      double rp2 = (double)(p2)*(double)(iy)/dratio / (double)n[1];
                      double theta_y = 2.0*PI*rp2;
                      std::complex<double> phase_y = std::complex<double>(cos(theta_y), sin(theta_y));
                      //if(p2 == n[1]/2) phase_y = std::real(phase_y);
                      //if(p2 == n[1]/2) phase_y = 1.0;

                      for(int izz = 0;izz < dimz_c;izz++) {
                          int p3 = offset_z + izz;
                          if(p3 > n[2]/2) p3 -= n[2];
                          //if(ct.sqrt_interpolation && (p3 == n[2]/2)) p3 = -p3;
                          double rp3 = (double)p3*(double)(iz)/dratio / (double)n[2];
                          double theta_z = 2.0*PI*rp3;
                          std::complex<double> phase_z = std::complex<double>(cos(theta_z), sin(theta_z));
                          //if(p3 == n[2]/2) phase_z = std::real(phase_z);
                          //if(p3 == n[2]/2) phase_z = 1.0;
                          std::complex<double> phase = phase_x * phase_y * phase_z;
                          shifted_coarse[idx] = phase * base_coarse[idx];
                          idx++;

                      }
                  }
              }

              // Backtransform phase shifted coefficients
              coarse_pwaves->FftInverse(shifted_coarse, backshifted_coarse);

              // Pack interpolated values into the fine grid
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      for(int izz = 0;izz < dimz_c;izz++) {
                          fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                              scale * backshifted_coarse[ixx*incx_c + iyy*incy_c + izz];
                      }
                  }
              }

          }

      }

  }

  // Loop over phase shifts.
  for(int ix = 0;ix < ratio;ix++) {

      for(int iy = 0;iy < ratio;iy++) {

          for(int iz = 0;iz < ratio;iz++) {

              // Multiply coarse transform by phase shifts and 
              int idx = 0;
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  int p1 = offset_x + ixx;
                  if(p1 > n[0]/2) p1 -= n[0];
                  if(ct.sqrt_interpolation && (p1 == n[0]/2)) p1 = -p1;
                  double rp1 = (double)(p1)*(double)(ix)/dratio / (double)n[0];
                  double theta_x = 2.0*PI*rp1;
                  std::complex<double> phase_x = std::complex<double>(cos(theta_x), sin(theta_x));
//                  if(p1 == n[0]/2) phase_x = std::real(phase_x);
                  //if(p1 == n[0]/2) phase_x = 1.0;

                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      int p2 = offset_y + iyy;
                      if(p2 > n[1]/2) p2 -= n[1];
                     if(ct.sqrt_interpolation && (p2 == n[1]/2)) p2 = -p2;
                      double rp2 = (double)(p2)*(double)(iy)/dratio / (double)n[1];
                      double theta_y = 2.0*PI*rp2;
                      std::complex<double> phase_y = std::complex<double>(cos(theta_y), sin(theta_y));
                      //if(p2 == n[1]/2) phase_y = std::real(phase_y);
                      //if(p2 == n[1]/2) phase_y = 1.0;

                      for(int izz = 0;izz < dimz_c;izz++) {
                          int p3 = offset_z + izz;
                          if(p3 > n[2]/2) p3 -= n[2];
                         if(ct.sqrt_interpolation && (p3 == n[2]/2)) p3 = -p3;
                          double rp3 = (double)p3*(double)(iz)/dratio / (double)n[2];
                          double theta_z = 2.0*PI*rp3;
                          std::complex<double> phase_z = std::complex<double>(cos(theta_z), sin(theta_z));
                          //if(p3 == n[2]/2) phase_z = std::real(phase_z);
                          //if(p3 == n[2]/2) phase_z = 1.0;
                          std::complex<double> phase = phase_x * phase_y * phase_z;
                          shifted_coarse[idx] = phase * base_coarse[idx];
                          idx++;

                      }
                  }
              }

              // Backtransform phase shifted coefficients
              coarse_pwaves->FftInverse(shifted_coarse, backshifted_coarse);

              // Pack interpolated values into the fine grid
              for(int ixx = 0;ixx < dimx_c;ixx++) {
                  for(int iyy = 0;iyy < dimy_c;iyy++) {
                      for(int izz = 0;izz < dimz_c;izz++) {
                          fine[(ratio*ixx + ix)*incx_f + (ratio*iyy + iy)*incy_f + (ratio*izz + iz)] = 
                              scale*backshifted_coarse[ixx*incx_c + iyy*incy_c + izz];

                      }
                  }
              }

          }

      }

  }


  delete [] backshifted_coarse;
  delete [] shifted_coarse;
  delete [] base_coarse;

}
