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
#include "main.h"
#include <math.h>
#include "transition.h"

Pw::Pw (BaseGrid &G, Lattice &L, int ratio, bool gamma_flag)
{

  // Grid parameters
  this->Grid = &G;
  this->comm = G.comm;
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
  this->distributed_plan = NULL;

  // Magnitudes of the g-vectors
  this->gmags = new double[this->pbasis]();

  // Mask array which is set to true for g-vectors with frequencies below the cutoff
  // and false otherwise.
  this->gmask = new bool[this->pbasis]();
  for(int i = 0;i < this->pbasis;i++) this->gmask[i] = false;

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
  int ivx = ratio * ((G.get_NX_GRID(1) - 1) / 2);
  int ivy = ratio * ((G.get_NY_GRID(1) - 1) / 2);
  int ivz = ratio * ((G.get_NZ_GRID(1) - 1) / 2);
  
  gvec[0] = (double)ivx * L.b0[0] + (double)ivy * L.b1[0] + (double)ivz * L.b2[0];
  gvec[0] *= L.celldm[0];
  gvec[1] = (double)ivx * L.b0[1] + (double)ivy * L.b1[1] + (double)ivz * L.b2[1];
  gvec[1] *= L.celldm[0];
  gvec[2] = (double)ivx * L.b0[2] + (double)ivy * L.b1[2] + (double)ivz * L.b2[2];
  gvec[2] *= L.celldm[0];

  this->gmax = 0.0;
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

              this->gmags[idx] = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
              this->gmax = std::max(this->gmax, this->gmags[idx]);
              g[idx].a[0] = gvec[0];
              g[idx].a[1] = gvec[1];
              g[idx].a[2] = gvec[2];

              if(abs(ivec[0]) > ivx) continue;
              if(abs(ivec[1]) > ivy) continue;
              if(abs(ivec[2]) > ivz) continue;

              // Gamma only exclude volume with x<0
              if((ivec[0] < 0) && this->is_gamma) continue;
              // Gamma only exclude plane with x = 0, y < 0
              if((ivec[0] == 0) && (ivec[1] < 0) && this->is_gamma) continue;
              // Gamma only exclude line with x = 0, y = 0, z < 0
              if((ivec[0] == 0) && (ivec[1] == 0) && (ivec[2] < 0) && this->is_gamma) continue;
              if((abs(ivec[0]) == this->global_dimx/2) || 
                 (abs(ivec[1]) == this->global_dimy/2) || 
                 (abs(ivec[2]) == this->global_dimz/2)) continue;

              this->gmask[idx] = true;
              this->ng++;
              //if((iy==0)&&(iz==0))printf("DDD  %d  %d  %d  %f\n",ivx,ix,ivec[0],gvec[0]);
          }
      }
  }

  MPI_Allreduce(MPI_IN_PLACE, &this->gmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  this->gcut = this->gmax;
  for(int idx = 0;idx < this->dimx*this->dimy*this->dimz;idx++)
  {
      if(this->gmags[idx] > this->gcut) {
          if(this->gmask[idx]) this->ng--;
          this->gmask[idx] = false;
      }
  }

  //int gcount = this->ng;
  //MPI_Allreduce(MPI_IN_PLACE, &gcount, 1, MPI_INT, MPI_SUM, comm);
  //printf("G-vector count  = %d\n", gcount);
  //printf("G-vector cutoff = %8.2f\n", sqrt(this->gcut));

  // Now set up plans
  if(G.get_NPES() == 1)
  {
      // Local cpu based fft plan(s). We use the array execute functions so the in and out arrays
      // here are dummies to enable the use of FFTW_MEASURE. The caller has to ensure alignment
      std::complex<double> *in = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * 2 * this->global_dimx*this->global_dimy*this->global_dimz);
      std::complex<double> *out = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * 2 * this->global_dimx*this->global_dimy*this->global_dimz);

      fftw_forward_plan = fftw_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), 
                FFTW_FORWARD, FFTW_MEASURE);

      fftw_backward_plan = fftw_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), 
                FFTW_BACKWARD, FFTW_MEASURE);

      fftw_forward_plan_inplace = fftw_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(in), 
                FFTW_FORWARD, FFTW_MEASURE);

      fftw_backward_plan_inplace = fftw_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(in), 
                FFTW_BACKWARD, FFTW_MEASURE);

      delete [] out;
      delete [] in;
  }
  else
  {
      // NPES > 1 so setup distributed fft plan
      int grid[3];

      // See if we can use pfft without remapping. In and out arrays must be equal in size.
      grid[0] = G.get_NX_GRID(ratio);
      grid[1] = G.get_NY_GRID(ratio);
      grid[2] = G.get_NZ_GRID(ratio);
      int dimx = G.get_PX0_GRID(ratio);
      int dimy = G.get_PY0_GRID(ratio);
      int dimz = G.get_PZ0_GRID(ratio);

      int pxoffset, pyoffset, pzoffset, nbuf, scaled=false, permute=0, usecollective=true;
      G.find_node_offsets(G.get_rank(), grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);
      distributed_plan = fft_3d_create_plan(comm,
                           grid[2], grid[1], grid[0],
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           scaled, permute, &nbuf, usecollective);

  }

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

int Pw::count_filtered_gvectors(double filter_factor)
{
  int pbasis = this->dimx * this->dimy * this->dimz;
  int gcount = 0;
  for(int idx=0;idx < pbasis;idx++)
  {
      if(this->gmags[idx] < filter_factor*this->gcut) gcount++;
  }
  MPI_Allreduce(MPI_IN_PLACE, &gcount, 1, MPI_INT, MPI_SUM, comm);
  return gcount;

}

void Pw::FftForward (double * in, std::complex<double> * out)
{
  std::complex<double> *buf = new std::complex<double>[pbasis];
  for(int i = 0;i < pbasis;i++) buf[i] = std::complex<double>(in[i], 0.0);

  if(Grid->get_NPES() == 1)
  {
      fftw_execute_dft (fftw_forward_plan,  reinterpret_cast<fftw_complex*>(buf), reinterpret_cast<fftw_complex*>(out));
  }
  else
  {
      fft_3d((FFT_DATA *)buf, (FFT_DATA *)out, -1, distributed_plan);
  }
  delete [] buf;
}


void Pw::FftForward (std::complex<double> * in, std::complex<double> * out)
{
    if(Grid->get_NPES() == 1)
    {
        if(in == out)
            fftw_execute_dft (fftw_forward_plan_inplace,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
        else
            fftw_execute_dft (fftw_forward_plan,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
    }
    else
    {
        if(in != out)
        {
            std::complex<double> *buf = new std::complex<double>[pbasis];
            for(int i = 0;i < pbasis;i++) buf[i] = in[i];
            fft_3d((FFT_DATA *)buf, (FFT_DATA *)out, -1, distributed_plan);
            delete [] buf;
        }
        else
        {
            fft_3d((FFT_DATA *)in, (FFT_DATA *)out, -1, distributed_plan);
        }
    }
}


void Pw::FftInverse (std::complex<double> * in, std::complex<double> * out)
{
    if(Grid->get_NPES() == 1)
    {
        if(in == out)
            fftw_execute_dft (fftw_backward_plan_inplace,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
        else
            fftw_execute_dft (fftw_backward_plan,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
    }
    else
    {
        if(in != out)
        {
            std::complex<double> *buf = new std::complex<double>[pbasis];
            for(int i = 0;i < pbasis;i++) buf[i] = in[i];
            fft_3d((FFT_DATA *)buf, (FFT_DATA *)out, 1, distributed_plan);
            delete [] buf;
        }
        else
        {
            fft_3d((FFT_DATA *)in, (FFT_DATA *)out, 1, distributed_plan);
        }
    }
}


Pw::~Pw(void)
{
  delete [] g;
  delete [] gmask;
  delete [] gmags;
  if(this->distributed_plan) fft_3d_destroy_plan(this->distributed_plan);
  if(Grid->get_NPES() == 1)
  {
      fftw_destroy_plan(fftw_backward_plan_inplace);
      fftw_destroy_plan(fftw_forward_plan_inplace);
      fftw_destroy_plan(fftw_backward_plan);
      fftw_destroy_plan(fftw_forward_plan);
  }
}

