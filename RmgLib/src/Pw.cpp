#if USE_PFFT
/*
 *
 * Copyright (c) 2015, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * @file
 *
 *
 * @section DESCRIPTION
 * Used to handle the integration of plane wave functionality into RMG.
 */

#include "Pw.h"
#include <math.h>

Pw::Pw (BaseGrid &G, Lattice &L, int ratio)
{

  // Grid parameters
  this->Grid = &G;
  this->L = &L;
  this->pbasis = G.get_P0_BASIS(ratio);
  this->hxgrid = G.get_hxgrid(ratio);
  this->hygrid = G.get_hygrid(ratio);
  this->hzgrid = G.get_hzgrid(ratio);
  this->dimx = G.get_PX0_GRID(ratio);
  this->dimy = G.get_PY0_GRID(ratio);
  this->dimz = G.get_PZ0_GRID(ratio);
  this->global_dimx = G.get_NX_GRID(ratio);
  this->global_dimy = G.get_NY_GRID(ratio);
  this->global_dimz = G.get_NZ_GRID(ratio);
  this->global_basis = this->global_dimx * this->global_dimy * this->global_dimz;

  this->gmags = new double[this->pbasis];
  this->gmask = new double[this->pbasis]();
  this->g = new gvector[this->pbasis];

  int idx=0;
  int ivec[3];
  double gvec[3];

  // Get cutoff
  ivec[0] = (this->dimx - 1) / 2;
  ivec[1] = (this->dimy - 1) / 2;
  ivec[2] = (this->dimz - 1) / 2;
  index_to_gvector(ivec, gvec);
  this->gcut = sqrt(gvec[0] * gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]) / 2.0;

  for(int ix = 0;ix < this->dimx;ix++) {
      for(int iy = 0;iy < this->dimy;iy++) {
          for(int iz = 0;iz < this->dimz;iz++) {
              ivec[0] = ix;
              ivec[1] = iy;
              ivec[2] = iz;
              index_to_gvector(ivec, gvec);
              this->gmags[idx] = sqrt(gvec[0] * gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]);
              if(this->gmags[idx] <= this->gcut) this->gmask[idx] = 1.0;
              this->g->a[0] = gvec[0];
              this->g->a[1] = gvec[1];
              this->g->a[2] = gvec[2];
              idx++;
          }
      }
  }

}

// Converts local index into fft array into corresponding unnormalized g-vector on global grid
void Pw::index_to_gvector(int *index, double *gvector)
{

  int ivector[3];

  this->Grid->find_node_offsets(this->Grid->get_rank(), this->global_dimx, this->global_dimy, this->global_dimz, &ivector[0], &ivector[1], &ivector[2]);

  // Compute fft permutations
  ivector[0] = ivector[0] + index[0];
  if(ivector[0] > this->global_dimx/2) ivector[0] = ivector[0] - this->global_dimx;
  if(ivector[0] == this->global_dimx/2) ivector[0] = 0;

  ivector[1] = ivector[1] + index[1];
  if(ivector[1] > this->global_dimy/2) ivector[1] = ivector[1] - this->global_dimy;
  if(ivector[1] == this->global_dimy/2) ivector[1] = 0;

  ivector[2] = ivector[2] + index[2];
  if(ivector[2] > this->global_dimz/2) ivector[2] = ivector[2] - this->global_dimz;
  if(ivector[2] == this->global_dimz/2) ivector[2] = 0;

  // Generate g-vectors from reciprocal lattice vectors
  gvector[0] = (double)ivector[0] * L->b0[0] + (double)ivector[1] * L->b1[0] + (double)ivector[2] * L->b2[0];
  gvector[1] = (double)ivector[0] * L->b0[1] + (double)ivector[1] * L->b1[1] + (double)ivector[2] * L->b2[1];
  gvector[2] = (double)ivector[0] * L->b0[2] + (double)ivector[1] * L->b1[2] + (double)ivector[2] * L->b2[2];
}


Pw::~Pw(void)
{
  delete [] g;
  delete [] gmask;
  delete [] gmags;
}

#endif
