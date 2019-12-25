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
#include <omp.h>
#include "Pw.h"
#include "transition.h"
#include "rmg_error.h"
#include "ErrorFuncs.h"

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
  this->distributed_plan.resize(ct.MG_THREADS_PER_NODE);
  this->distributed_plan_f.resize(ct.MG_THREADS_PER_NODE);
  for(int thread=0;thread < ct.MG_THREADS_PER_NODE;thread++) this->distributed_plan[thread] = NULL;
  for(int thread=0;thread < ct.MG_THREADS_PER_NODE;thread++) this->distributed_plan_f[thread] = NULL;

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

  //            if(abs(ivec[0]) > ivx) continue;
  //            if(abs(ivec[1]) > ivy) continue;
  //            if(abs(ivec[2]) > ivz) continue;

              // Gamma only exclude volume with x<0
              if((ivec[0] < 0) && this->is_gamma) continue;
              // Gamma only exclude plane with x = 0, y < 0
              if((ivec[0] == 0) && (ivec[1] < 0) && this->is_gamma) continue;
              // Gamma only exclude line with x = 0, y = 0, z < 0
              if((ivec[0] == 0) && (ivec[1] == 0) && (ivec[2] < 0) && this->is_gamma) continue;
   //           if((abs(ivec[0]) == this->global_dimx/2) || 
    //             (abs(ivec[1]) == this->global_dimy/2) || 
    //             (abs(ivec[2]) == this->global_dimz/2)) continue;

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
      std::complex<double> *in = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * this->global_dimx*this->global_dimy*this->global_dimz);
      std::complex<double> *out = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * this->global_dimx*this->global_dimy*this->global_dimz);

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

      fftwf_forward_plan = fftwf_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out), 
                FFTW_FORWARD, FFTW_MEASURE);

      fftwf_backward_plan = fftwf_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out), 
                FFTW_BACKWARD, FFTW_MEASURE);

      fftwf_forward_plan_inplace = fftwf_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(in), 
                FFTW_FORWARD, FFTW_MEASURE);

      fftwf_backward_plan_inplace = fftwf_plan_dft_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(in), 
                FFTW_BACKWARD, FFTW_MEASURE);

      delete [] out;
      delete [] in;

#if GPU_ENABLED
      num_streams = ct.OMP_THREADS_PER_NODE;
      // Gpu streams and plans
      streams.resize(num_streams);
      gpu_plans.resize(num_streams);
      gpu_plans_f.resize(num_streams);
      host_bufs.resize(num_streams);
      dev_bufs.resize(num_streams);

      for (int i = 0; i < num_streams; i++)
          RmgCudaError(__FILE__, __LINE__, cudaStreamCreateWithFlags(&streams[i],cudaStreamNonBlocking), "Problem creating cuda stream.");

      for (int i = 0; i < num_streams; i++)
      {
         cufftPlan3d(&gpu_plans[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_Z2Z);
         cufftPlan3d(&gpu_plans_f[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_C2C);
         cufftSetStream(gpu_plans[i], streams[i]);
         cufftSetStream(gpu_plans_f[i], streams[i]);
         RmgCudaError(__FILE__, __LINE__, 
             cudaMallocHost((void **)&host_bufs[i],  this->global_dimx*this->global_dimy*this->global_dimz * sizeof(std::complex<double>)),
             "Error: cudaMallocHost failed.\n");

         RmgCudaError(__FILE__, __LINE__, 
             cudaMalloc((void **)&dev_bufs[i],  this->global_dimx*this->global_dimy*this->global_dimz * sizeof(std::complex<double>)),
             "Error: cudaMalloc failed.\n");
      }
 
#endif

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

      int pxoffset, pyoffset, pzoffset, nbuf, scaled=false, permute=0, usecollective=false;
      G.find_node_offsets(G.get_rank(), grid[0], grid[1], grid[2], &pxoffset, &pyoffset, &pzoffset);
      for(int thread=0;thread < ct.MG_THREADS_PER_NODE;thread++)
      {
          distributed_plan[thread] = fft_3d_create_plan<fftw_complex, double>(comm,
                           grid[2], grid[1], grid[0],
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           pzoffset, pzoffset + dimz - 1,
                           pyoffset, pyoffset + dimy - 1,
                           pxoffset, pxoffset + dimx - 1,
                           scaled, permute, &nbuf, usecollective);

          distributed_plan_f[thread] = fft_3d_create_plan<fftwf_complex, float>(comm,
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
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if GPU_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      for(int i = 0;i < pbasis;i++) tptr[i] = std::complex<double>(in[i], 0.0);
      cudaMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyHostToDevice, streams[tid]);
      cufftExecZ2Z(gpu_plans[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid], CUFFT_FORWARD);
      cudaMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyDeviceToHost, streams[tid]);
      cudaStreamSynchronize(streams[tid]);
      for(int i = 0;i < pbasis;i++) out[i] = tptr[i];
#else
      std::complex<double> *buf = new std::complex<double>[pbasis];
      for(int i = 0;i < pbasis;i++) buf[i] = std::complex<double>(in[i], 0.0);
      fftw_execute_dft (fftw_forward_plan,  reinterpret_cast<fftw_complex*>(buf), reinterpret_cast<fftw_complex*>(out));
      delete [] buf;
#endif
  }
  else
  {
      std::complex<double> *buf = new std::complex<double>[pbasis];
      for(int i = 0;i < pbasis;i++) buf[i] = std::complex<double>(in[i], 0.0);
      fft_3d((fftw_complex *)buf, (fftw_complex *)out, -1, distributed_plan[tid]);
      delete [] buf;
  }
}

void Pw::FftForward (float * in, std::complex<float> * out)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  std::complex<float> *buf = new std::complex<float>[pbasis];
  for(int i = 0;i < pbasis;i++) buf[i] = std::complex<float>(in[i], 0.0);
  
  if(Grid->get_NPES() == 1)
  {   
      fftwf_execute_dft (fftwf_forward_plan,  reinterpret_cast<fftwf_complex*>(buf), reinterpret_cast<fftwf_complex*>(out));
  }
  else
  {   
      fft_3d((fftwf_complex *)buf, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
  }
  delete [] buf;
}

void Pw::FftForward (std::complex<double> * in, std::complex<double> * out)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if GPU_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      for(int i = 0;i < pbasis;i++) tptr[i] = in[i];
      cudaMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyHostToDevice, streams[tid]);
      cufftExecZ2Z(gpu_plans[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid], CUFFT_FORWARD);
      cudaMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyDeviceToHost, streams[tid]);
      cudaStreamSynchronize(streams[tid]);
      for(int i = 0;i < pbasis;i++) out[i] = tptr[i];
#else
      if(in == out)
          fftw_execute_dft (fftw_forward_plan_inplace,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
      else
          fftw_execute_dft (fftw_forward_plan,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
#endif
  }
  else
  {
      if(in != out)
      {
          std::complex<double> *buf = new std::complex<double>[pbasis];
          for(int i = 0;i < pbasis;i++) buf[i] = in[i];
          fft_3d((fftw_complex *)buf, (fftw_complex *)out, -1, distributed_plan[tid]);
          delete [] buf;
      }
      else
      {
          fft_3d((fftw_complex *)in, (fftw_complex *)out, -1, distributed_plan[tid]);
      }
  }
}


void Pw::FftForward (std::complex<float> * in, std::complex<float> * out)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  if(Grid->get_NPES() == 1)
  {   
#if GPU_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      for(int i = 0;i < pbasis;i++) tptr[i] = in[i];
      cudaMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<float>), cudaMemcpyHostToDevice, streams[tid]);
      cufftExecC2C(gpu_plans_f[tid], (cufftComplex*)dev_bufs[tid], (cufftComplex*)dev_bufs[tid], CUFFT_FORWARD);
      cudaMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<float>), cudaMemcpyDeviceToHost, streams[tid]);
      cudaStreamSynchronize(streams[tid]);
      for(int i = 0;i < pbasis;i++) out[i] = tptr[i];
#else
      if(in == out)
          fftwf_execute_dft (fftwf_forward_plan_inplace,  reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
      else
          fftwf_execute_dft (fftwf_forward_plan,  reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
#endif
  }
  else
  {   
      if(in != out)
      {   
          std::complex<float> *buf = new std::complex<float>[pbasis];
          for(int i = 0;i < pbasis;i++) buf[i] = in[i];
          fft_3d((fftwf_complex *)buf, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
          delete [] buf;
      }
      else
      {
          fft_3d((fftwf_complex *)in, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
      }
  }
}


void Pw::FftInverse (std::complex<double> * in, std::complex<double> * out)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if GPU_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      for(int i = 0;i < pbasis;i++) tptr[i] = in[i];
      cudaMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyHostToDevice, streams[tid]);
      cufftExecZ2Z(gpu_plans[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid], CUFFT_INVERSE);
      cudaMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<double>), cudaMemcpyDeviceToHost, streams[tid]);
      cudaStreamSynchronize(streams[tid]);
      for(int i = 0;i < pbasis;i++) out[i] = tptr[i];
#else
      if(in == out)
          fftw_execute_dft (fftw_backward_plan_inplace,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
      else
          fftw_execute_dft (fftw_backward_plan,  reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
#endif
  }
  else
  {
      if(in != out)
      {
          std::complex<double> *buf = new std::complex<double>[pbasis];
          for(int i = 0;i < pbasis;i++) buf[i] = in[i];
          fft_3d((fftw_complex *)buf, (fftw_complex *)out, 1, distributed_plan[tid]);
          delete [] buf;
      }
      else
      {
          fft_3d((fftw_complex *)in, (fftw_complex *)out, 1, distributed_plan[tid]);
      }
  }
}

void Pw::FftInverse (std::complex<float> * in, std::complex<float> * out)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  if(Grid->get_NPES() == 1)
  {
#if GPU_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      for(int i = 0;i < pbasis;i++) tptr[i] = in[i];
      cudaMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<float>), cudaMemcpyHostToDevice, streams[tid]);
      cufftExecC2C(gpu_plans_f[tid], (cufftComplex*)dev_bufs[tid], (cufftComplex*)dev_bufs[tid], CUFFT_INVERSE);
      cudaMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<float>), cudaMemcpyDeviceToHost, streams[tid]);
      cudaStreamSynchronize(streams[tid]);
      for(int i = 0;i < pbasis;i++) out[i] = tptr[i];
#else 
      if(in == out)
      {
          fftwf_execute_dft (fftwf_backward_plan_inplace,  reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
      }
      else
      {
          std::complex<float> *buf = new std::complex<float>[pbasis];
          for(int i = 0;i < pbasis;i++) buf[i] = in[i];
          fftwf_execute_dft (fftwf_backward_plan,  (fftwf_complex *)buf, reinterpret_cast<fftwf_complex*>(out));
          delete [] buf;
      }
#endif
  }
  else
  {   
      if(in != out)
      {   
          std::complex<float> *buf = new std::complex<float>[pbasis];
          for(int i = 0;i < pbasis;i++) buf[i] = in[i];
          fft_3d((fftwf_complex *)buf, (fftwf_complex *)out, 1, distributed_plan_f[tid]);
          delete [] buf;
      }
      else
      {   
          fft_3d((fftwf_complex *)in, (fftwf_complex *)out, 1, distributed_plan_f[tid]);
      }
  }
}

Pw::~Pw(void)
{
  delete [] g;
  delete [] gmask;
  delete [] gmags;

  for(int thread=0;thread < ct.MG_THREADS_PER_NODE;thread++)
  {
      if(this->distributed_plan[thread]) fft_3d_destroy_plan(this->distributed_plan[thread]);
      if(this->distributed_plan_f[thread]) fft_3d_destroy_plan(this->distributed_plan_f[thread]);
  }

  if(Grid->get_NPES() == 1)
  {
      fftw_destroy_plan(fftw_backward_plan_inplace);
      fftw_destroy_plan(fftw_forward_plan_inplace);
      fftw_destroy_plan(fftw_backward_plan);
      fftw_destroy_plan(fftw_forward_plan);
      fftwf_destroy_plan(fftwf_backward_plan_inplace);
      fftwf_destroy_plan(fftwf_forward_plan_inplace);
      fftwf_destroy_plan(fftwf_backward_plan);
      fftwf_destroy_plan(fftwf_forward_plan);
#if GPU_ENABLED
      for (int i = 0; i < num_streams; i++)
      {
         cufftDestroy(gpu_plans[i]);
         cufftDestroy(gpu_plans_f[i]);
         RmgCudaError(__FILE__, __LINE__, cudaFreeHost(host_bufs[i]), "Error: cudaFreeHost failed.\n");
         RmgCudaError(__FILE__, __LINE__, cudaFree(dev_bufs[i]), "Error: cudaFreeHost failed.\n");
      }

      // Gpu streams and plans
      for (int i = 0; i < num_streams; i++)
          RmgCudaError(__FILE__, __LINE__, cudaStreamDestroy(streams[i]), "Problem freeing cuda stream.");


#endif
  }
}

