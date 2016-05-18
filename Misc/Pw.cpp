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

#include "Pw.h"
#include <math.h>

Pw::Pw (BaseGrid &G, Lattice &L, int ratio, bool gamma_flag, MPI_Comm comm)
{

  // Grid parameters
  this->Grid = &G;
  this->L = &L;
  this->is_gamma = gamma_flag;
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

  // Magnitudes of the g-vectors
  this->gmags = new double[this->pbasis]();

  // Mask array which is set to 1.0 for g-vectors with frequencies below the cutoff
  // and 0.0 otherwise.
  this->gmask = new double[this->pbasis]();

  // G-vector storage. Vectors with frequencies above the cutoff are stored as zeros
  this->g = new gvector[this->pbasis];

  // Number of g-vectors
  this->ng = 0;

  int ivec[3];
  double gvec[3];

  // zero out g-vector array
  for(int ix=0;ix < this->pbasis;ix++) {
      this->g->a[0] = 0.0;
      this->g->a[1] = 0.0;
      this->g->a[2] = 0.0;
  }

  // Get G^2 cutoff.
  ivec[0] = (this->global_dimx - 1) / 2;
  ivec[1] = (this->global_dimy - 1) / 2;
  ivec[2] = (this->global_dimz - 1) / 2;
  gvec[0] = (double)ivec[0] * L.b0[0] + (double)ivec[1] * L.b1[0] + (double)ivec[2] * L.b2[0];
  gvec[0] *= L.celldm[0];
  gvec[1] = (double)ivec[0] * L.b0[1] + (double)ivec[1] * L.b1[1] + (double)ivec[2] * L.b2[1];
  gvec[1] *= L.celldm[0];
  gvec[2] = (double)ivec[0] * L.b0[2] + (double)ivec[1] * L.b1[2] + (double)ivec[2] * L.b2[2];
  gvec[2] *= L.celldm[0];

  this->gmax = gvec[0] * gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
  this->gcut = this->gmax / 3.0;
  
  int idx = -1;
  for(int ix = 0;ix < this->dimx;ix++) {
      for(int iy = 0;iy < this->dimy;iy++) {
          for(int iz = 0;iz < this->dimz;iz++) {
              idx++;
              ivec[0] = ix;
              ivec[1] = iy;
              ivec[2] = iz;

              // On input local index is stored in ivec. On output global index
              // is stored in ivec and the associated g-vector is stored in g
              index_to_gvector(ivec, gvec);
              // Gamma only exclude volume with x<0
              if((ivec[0] < 0) && this->is_gamma) continue;
              // Gamma only exclude plane with x = 0, y < 0
              if((ivec[0] == 0) && (ivec[1] < 0) && this->is_gamma) continue;
              // Gamma only exclude line with x = 0, y = 0, z < 0
              if((ivec[0] == 0) && (ivec[1] == 0) && (ivec[2] < 0) && this->is_gamma) continue;

              this->gmags[idx] = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
              g[idx].a[0] = gvec[0];
              g[idx].a[1] = gvec[1];
              g[idx].a[2] = gvec[2];

              if(((int)fabs(round(gvec[0])) == this->global_dimx/2) || 
                 ((int)fabs(round(gvec[1])) == this->global_dimy/2) || 
                 ((int)fabs(round(gvec[2])) == this->global_dimz/2)) continue;

              if(this->gmags[idx] <= this->gcut) {
                  this->gmask[idx] = 1.0;
                  this->ng++;
              }
          }
      }
  }

  //printf("G-vector count  = %d\n", this->ng);
  //printf("G-vector cutoff = %8.2f\n", this->gcut);

}

// Converts local index into fft array into corresponding g-vector on global grid
// On input index holds the local index
// On output index holds the global index and gvector holds the corresponding g-vector
void Pw::index_to_gvector(int *index, double *gvector)
{

  int ivector[3];

  this->Grid->find_node_offsets(this->Grid->get_rank(), this->global_dimx, this->global_dimy, this->global_dimz, &ivector[0], &ivector[1], &ivector[2]);

  // Compute fft permutations
  ivector[0] = ivector[0] + index[0];
  if(ivector[0] > this->global_dimx/2) ivector[0] = ivector[0] - this->global_dimx;

  ivector[1] = ivector[1] + index[1];
  if(ivector[1] > this->global_dimy/2) ivector[1] = ivector[1] - this->global_dimy;

  ivector[2] = ivector[2] + index[2];
  if(ivector[2] > this->global_dimz/2) ivector[2] = ivector[2] - this->global_dimz;

  // Generate g-vectors from reciprocal lattice vectors
  gvector[0] = (double)ivector[0] * L->b0[0] + (double)ivector[1] * L->b1[0] + (double)ivector[2] * L->b2[0];
  gvector[1] = (double)ivector[0] * L->b0[1] + (double)ivector[1] * L->b1[1] + (double)ivector[2] * L->b2[1];
  gvector[2] = (double)ivector[0] * L->b0[2] + (double)ivector[1] * L->b1[2] + (double)ivector[2] * L->b2[2];

  gvector[0] *= L->celldm[0];
  gvector[1] *= L->celldm[0];
  gvector[2] *= L->celldm[0];

  index[0] = ivector[0];
  index[1] = ivector[1];
  index[2] = ivector[2];
}


Pw::~Pw(void)
{
  delete [] g;
  delete [] gmask;
  delete [] gmags;
}

