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
#include "GpuAlloc.h"

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
  this->global_basis = (size_t)this->global_dimx * (size_t)this->global_dimy * (size_t)this->global_dimz;
  this->global_basis_alloc = (size_t)this->global_dimx * (size_t)this->global_dimy * (size_t)(2*(this->global_dimz/2 + 1));
  this->global_basis_packed = (size_t)this->global_dimx * (size_t)this->global_dimy * (size_t)(this->global_dimz/2 + 1);
  int max_threads = std::max(ct.MG_THREADS_PER_NODE, ct.OMP_THREADS_PER_NODE);
  this->distributed_plan.resize(max_threads);
  this->distributed_plan_f.resize(max_threads);
  for(int thread=0;thread < max_threads;thread++) this->distributed_plan[thread] = NULL;
  for(int thread=0;thread < max_threads;thread++) this->distributed_plan_f[thread] = NULL;

  // Magnitudes of the g-vectors
  this->gmags = new double[this->pbasis]();

  // Mask array which is set to true for g-vectors with frequencies below the cutoff
  // and false otherwise.
  this->gmask = new bool[this->pbasis]();
  for(size_t i = 0;i < this->pbasis;i++) this->gmask[i] = false;

  // G-vector storage. Vectors with frequencies above the cutoff are stored as zeros
  this->g = new gvector[this->pbasis];

  // Number of g-vectors
  this->ng = 0;

  int ivec[3];
  double gvec[3];

  // zero out g-vector array
  for(size_t ix=0;ix < this->pbasis;ix++) {
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
  int64_t idx = -1;
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

  if(L.get_ibrav_type() == CUBIC_FC) gmax /= sqrt(2.0);
  if(L.get_ibrav_type() == CUBIC_BC) gmax /= sqrt(6.0);

  this->gcut = this->gmax;
  for(size_t idx = 0;idx < this->pbasis;idx++)
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

#if CUDA_ENABLED || HIP_ENABLED
      num_streams = ct.OMP_THREADS_PER_NODE;
      num_streams = std::max(ct.MG_THREADS_PER_NODE, ct.OMP_THREADS_PER_NODE);
      // Gpu streams and plans
      streams.resize(num_streams);
      gpu_plans.resize(num_streams);
      gpu_plans_f.resize(num_streams);
      gpu_plans_r2c.resize(num_streams);
      gpu_plans_c2r.resize(num_streams);
      gpu_plans_d2z.resize(num_streams);
      gpu_plans_z2d.resize(num_streams);
      host_bufs.resize(num_streams);
      dev_bufs.resize(num_streams);

#endif

#if CUDA_ENABLED
      for (int i = 0; i < num_streams; i++)
          RmgGpuError(__FILE__, __LINE__, gpuStreamCreateWithFlags(&streams[i], gpuStreamNonBlocking), "Problem creating gpu stream.");

      for (int i = 0; i < num_streams; i++)
      {
         cufftResult res;
         res = cufftPlan3d(&gpu_plans[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_Z2Z);
         if(res != CUFFT_SUCCESS) printf("ERROR making Z2Z plan\n");

         res = cufftPlan3d(&gpu_plans_f[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_C2C);
         if(res != CUFFT_SUCCESS) printf("ERROR making C2C plan\n");

         res = cufftPlan3d(&gpu_plans_d2z[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_D2Z);
         if(res != CUFFT_SUCCESS) printf("ERROR making D2Z plan\n");

         res = cufftPlan3d(&gpu_plans_z2d[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_Z2D);
         if(res != CUFFT_SUCCESS) printf("ERROR making Z2D plan\n");

         res = cufftPlan3d(&gpu_plans_r2c[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_R2C);
         if(res != CUFFT_SUCCESS) printf("ERROR making R2C plan\n");

         res = cufftPlan3d(&gpu_plans_c2r[i], this->global_dimx, this->global_dimy, this->global_dimz, CUFFT_C2R);
         if(res != CUFFT_SUCCESS) printf("ERROR making C2R plan\n");

         cufftSetStream(gpu_plans[i], streams[i]);
         cufftSetStream(gpu_plans_f[i], streams[i]);
         cufftSetStream(gpu_plans_d2z[i], streams[i]);
         cufftSetStream(gpu_plans_z2d[i], streams[i]);
         cufftSetStream(gpu_plans_r2c[i], streams[i]);
         cufftSetStream(gpu_plans_c2r[i], streams[i]);

         RmgGpuError(__FILE__, __LINE__, 
             gpuMallocHost((void **)&host_bufs[i],  this->global_basis_alloc * sizeof(std::complex<double>)),
             "Error: gpuMallocHost failed.\n");

         RmgGpuError(__FILE__, __LINE__, 
             gpuMalloc((void **)&dev_bufs[i],  this->global_basis_alloc * sizeof(std::complex<double>)),
             "Error: gpuMalloc failed.\n");
      }
#elif HIP_ENABLED
      if(num_streams > 1)
      {
          for (int i = 0; i < num_streams; i++)
              RmgGpuError(__FILE__, __LINE__, gpuStreamCreateWithFlags(&streams[i], gpuStreamNonBlocking), "Problem creating gpu stream.");
      }
      else
      {
          // Only 1 thread so use default stream
          streams[0] = 0;
      }

      work_bufs.resize(num_streams);
      roc_x_info.resize(num_streams);
      gpu_plans_inv.resize(num_streams);
      gpu_plans_f_inv.resize(num_streams);

      size_t lengths[3], work_size=0, tmp_size;
      rocfft_status status;
      lengths[0] = this->global_dimz;
      lengths[1] = this->global_dimy;
      lengths[2] = this->global_dimx;

      for (int i = 0; i < num_streams; i++)
      {
          status = rocfft_plan_create(&gpu_plans[i], rocfft_placement_inplace, 
                   rocfft_transform_type_complex_forward, rocfft_precision_double, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making Z2Z plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making Z2Z plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_inv[i], rocfft_placement_inplace, 
                   rocfft_transform_type_complex_inverse, rocfft_precision_double, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making Z2Z plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_inv[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making Z2Z plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_f[i], rocfft_placement_inplace, 
                   rocfft_transform_type_complex_forward, rocfft_precision_single, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making C2C plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_f[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making C2C plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_f_inv[i], rocfft_placement_inplace, 
                   rocfft_transform_type_complex_inverse, rocfft_precision_single, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making C2C plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_f_inv[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making C2C plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_d2z[i], rocfft_placement_notinplace, 
                   rocfft_transform_type_real_forward, rocfft_precision_double, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making D2Z plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_d2z[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making D2Z plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_z2d[i], rocfft_placement_notinplace, 
                   rocfft_transform_type_real_inverse, rocfft_precision_double, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making Z2D plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_z2d[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making Z2D plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_r2c[i], rocfft_placement_inplace, 
                   rocfft_transform_type_real_forward, rocfft_precision_single, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making R2C plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_r2c[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making R2C plans.\n");
          work_size = std::max(work_size, tmp_size);

          status = rocfft_plan_create(&gpu_plans_c2r[i], rocfft_placement_inplace, 
                   rocfft_transform_type_real_inverse, rocfft_precision_single, 3, lengths, 1, NULL);
          if(status != rocfft_status_success) printf("ERROR making C2R plans.\n");
          status = rocfft_plan_get_work_buffer_size(gpu_plans_c2r[i], &tmp_size);
          if(status != rocfft_status_success) printf("ERROR making C2R plans.\n");
          work_size = std::max(work_size, tmp_size);


          status = rocfft_execution_info_create(&roc_x_info[i]);
          if(status != rocfft_status_success) printf("ERROR making rocfft plans.\n");
          status = rocfft_execution_info_set_stream(roc_x_info[i], streams[i]);
          if(status != rocfft_status_success) printf("ERROR setting rocfft streams.\n");
          if(work_size)
          {
              RmgGpuError(__FILE__, __LINE__, 
                         gpuMalloc((void **)&work_bufs[i], work_size),
                         "Error: gpuMalloc failed.\n");
              status = rocfft_execution_info_set_work_buffer(roc_x_info[i], work_bufs[i], work_size);
          }


          RmgGpuError(__FILE__, __LINE__, 
                   gpuMallocHost((void **)&host_bufs[i],  this->global_basis_alloc * sizeof(std::complex<double>)),
                   "Error: gpuMallocHost failed.\n");

          RmgGpuError(__FILE__, __LINE__, 
                   gpuMalloc((void **)&dev_bufs[i],  this->global_basis_alloc * sizeof(std::complex<double>)),
                   "Error: gpuMalloc failed.\n");

      }

#else
      // Local cpu based fft plan(s). We use the array execute functions so the in and out arrays
      // here are dummies to enable the use of FFTW_MEASURE. The caller has to ensure alignment
      std::complex<double> *in = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * this->global_basis_alloc);
      std::complex<double> *out = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * this->global_basis_alloc);

      fftw_r2c_forward_plan = fftw_plan_dft_r2c_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     (double *)in, reinterpret_cast<fftw_complex*>(out), FFTW_MEASURE);

      fftw_r2c_backward_plan = fftw_plan_dft_c2r_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), (double *)out, FFTW_MEASURE);

      fftw_r2c_forward_plan_inplace = fftw_plan_dft_r2c_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     (double *)in, reinterpret_cast<fftw_complex*>(in), FFTW_MEASURE);

      fftw_r2c_backward_plan_inplace = fftw_plan_dft_c2r_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftw_complex*>(in), (double *)in, FFTW_MEASURE);

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


      fftwf_r2c_forward_plan = fftwf_plan_dft_r2c_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     (float *)in, reinterpret_cast<fftwf_complex*>(out), FFTW_PATIENT);

      fftwf_r2c_backward_plan = fftwf_plan_dft_c2r_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), (float *)out, FFTW_PATIENT);

      fftwf_r2c_forward_plan_inplace = fftwf_plan_dft_r2c_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     (float *)in, reinterpret_cast<fftwf_complex*>(in), FFTW_PATIENT);

      fftwf_r2c_backward_plan_inplace = fftwf_plan_dft_c2r_3d (this->global_dimx, this->global_dimy, this->global_dimz, 
                     reinterpret_cast<fftwf_complex*>(in), (float *)in, FFTW_PATIENT);


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

      fftw_free(out);
      fftw_free(in);
      int nthreads = std::max(ct.OMP_THREADS_PER_NODE, ct.MG_THREADS_PER_NODE);
      host_bufs.resize(nthreads);
      for (int i = 0; i < nthreads; i++)
      {
          host_bufs[i] = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * this->global_basis_alloc);
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

size_t Pw::count_filtered_gvectors(double filter_factor)
{
  size_t pbasis = (size_t)this->dimx * (size_t)this->dimy * (size_t)this->dimz;
  size_t gcount = 0;
  for(size_t idx=0;idx < pbasis;idx++)
  {
      if(this->gmags[idx] < filter_factor*this->gcut) gcount++;
  }
  MPI_Allreduce(MPI_IN_PLACE, &gcount, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
  return gcount;

}

void Pw::FftForward (double * in, std::complex<double> * out)
{
    FftForward(in, out, true, true, true);
}

void Pw::FftForward (double * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      double *tptr = (double *)host_bufs[tid];
      if(copy_to_dev)
      {
          if(tptr != in)  for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
	  gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(double), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecD2Z(gpu_plans_d2z[tid], (cufftDoubleReal*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid]);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
	  gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], global_basis_alloc*sizeof(std::complex<double>), gpuMemcpyDeviceToHost, streams[tid]);
	  gpuStreamSynchronize(streams[tid]);
          std::complex<double> *tptr1 = (std::complex<double> *)host_bufs[tid];
	  for(size_t i = 0;i < pbasis;i++) out[i] = tptr1[i];
      }
#elif HIP_ENABLED
      double *tptr = (double *)host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
	  if(tptr != in) for(size_t i = 0;i < global_basis_alloc;i++) tptr[i] = in[i];
	  hipMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(double), gpuMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_d2z[tid], (void**) &dev_bufs[tid], (void**) &dev_bufs[tid], roc_x_info[tid]);
      if(copy_from_dev)
      {
	  hipMemcpyAsync(host_bufs[tid], dev_bufs[tid], global_basis_alloc*sizeof(std::complex<double>), gpuMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          std::complex<double> *tptr1 = (std::complex<double> *)host_bufs[tid];
	  for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
#else
    if((void *)in == (void *)out)
    {
        fftw_execute_dft_r2c (fftw_r2c_forward_plan_inplace,  in, reinterpret_cast<fftw_complex*>(out));
    }
    else
    {
        fftw_execute_dft_r2c (fftw_r2c_forward_plan,  in, reinterpret_cast<fftw_complex*>(out));
    }
#endif
  }
  else
  {
      std::complex<double> *buf = new std::complex<double>[pbasis];
      for(size_t i = 0;i < pbasis;i++) buf[i] = std::complex<double>(in[i], 0.0);
      fft_3d((fftw_complex *)buf, (fftw_complex *)out, -1, distributed_plan[tid]);
      delete [] buf;
  }
}

void Pw::FftForward (float * in, std::complex<float> * out)
{
    FftForward(in, out, true, true, true);
}

void Pw::FftForward (float * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  
  if(Grid->get_NPES() == 1)
  {   
#if CUDA_ENABLED
      float *tptr = (float *)host_bufs[tid];
      if(copy_to_dev)
      {
	  if(tptr != in) for(size_t i = 0;i < global_basis_alloc;i++) tptr[i] = in[i];
          gpuStreamSynchronize(streams[tid]);
	  gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(float), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecR2C(gpu_plans_r2c[tid], (cufftReal*)dev_bufs[tid], (cufftComplex*)dev_bufs[tid]);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
	  gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], global_basis_alloc*sizeof(std::complex<float>), gpuMemcpyDeviceToHost, streams[tid]);
	  gpuStreamSynchronize(streams[tid]);
          std::complex<float> *tptr1 = (std::complex<float> *)host_bufs[tid];
	  for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
#elif HIP_ENABLED
      float *tptr = (float *)host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
	  if(tptr != in) for(size_t i = 0;i < global_basis_alloc;i++) tptr[i] = in[i];
	  hipMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(float), gpuMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_r2c[tid], (void**) &dev_bufs[tid], (void**) &dev_bufs[tid], roc_x_info[tid]);
      if(copy_from_dev)
      {
	  hipMemcpyAsync(host_bufs[tid], dev_bufs[tid], global_basis_alloc*sizeof(std::complex<float>), gpuMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          std::complex<float> *tptr1 = (std::complex<float> *)host_bufs[tid];
	  for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
#else
      if((void *)in == (void *)out)
      {
          fftwf_execute_dft_r2c (fftwf_r2c_forward_plan_inplace,  in, reinterpret_cast<fftwf_complex*>(out));
      }
      else
      {
          fftwf_execute_dft_r2c (fftwf_r2c_forward_plan,  in, reinterpret_cast<fftwf_complex*>(out));
      }
#endif
  }
  else
  {   
      std::complex<float> *buf = new std::complex<float>[pbasis];
      for(size_t i = 0;i < pbasis;i++) buf[i] = std::complex<float>(in[i], 0.0);
      fft_3d((fftwf_complex *)buf, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
      delete [] buf;
  }
}

void Pw::FftForward (std::complex<double> * in, std::complex<double> * out)
{
    FftForward(in, out, true, true, true);
}

void Pw::FftForward (std::complex<double> * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecZ2Z(gpu_plans[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid], CUFFT_FORWARD);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
          gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<double>), gpuMemcpyDeviceToHost, streams[tid]);
          gpuStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
#elif HIP_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], tptr, pbasis*sizeof(std::complex<double>), hipMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans[tid], (void**) &dev_bufs[tid], NULL, roc_x_info[tid]);
      if(copy_from_dev)
      {   
          hipMemcpyAsync(tptr, dev_bufs[tid], pbasis*sizeof(std::complex<double>), hipMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
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
          for(size_t i = 0;i < pbasis;i++) buf[i] = in[i];
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
    FftForward(in, out, true, true, true);
}

void Pw::FftForward (std::complex<float> * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  if(Grid->get_NPES() == 1)
  {   
#if CUDA_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<float>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecC2C(gpu_plans_f[tid], (cufftComplex*)dev_bufs[tid], (cufftComplex*)dev_bufs[tid], CUFFT_FORWARD);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
          gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<float>), gpuMemcpyDeviceToHost, streams[tid]);
          gpuStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
#elif HIP_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], tptr, pbasis*sizeof(std::complex<float>), hipMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_f[tid], (void**) &dev_bufs[tid], NULL, roc_x_info[tid]);
      if(copy_from_dev)
      {   
          hipMemcpyAsync(tptr, dev_bufs[tid], pbasis*sizeof(std::complex<float>), hipMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
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
          for(size_t i = 0;i < pbasis;i++) buf[i] = in[i];
          fft_3d((fftwf_complex *)buf, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
          delete [] buf;
      }
      else
      {
          fft_3d((fftwf_complex *)in, (fftwf_complex *)out, -1, distributed_plan_f[tid]);
      }
  }
}

void Pw::FftInverse (std::complex<double> * in, double * out)
{
    FftInverse(in, out, true, true, true);
}

void Pw::FftInverse (std::complex<double> * in, double * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{

  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      std::complex<double> *tptr = (std::complex<double> *)host_bufs[tid];
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecZ2D(gpu_plans_z2d[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleReal*)dev_bufs[tid]);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
          gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], global_basis_alloc*sizeof(double), gpuMemcpyDeviceToHost, streams[tid]);
          gpuStreamSynchronize(streams[tid]);
          double *tptr1 = (double *)host_bufs[tid];
          if(out != tptr1) for(size_t i = 0;i < pbasis;i++) out[i] = tptr1[i];
      }
#elif HIP_ENABLED
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          std::complex<double> *tptr = (std::complex<double> *)host_bufs[tid];
          for(size_t i = 0;i < global_basis_alloc;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(std::complex<double>), gpuMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_z2d[tid], (void**) &dev_bufs[tid], (void**) &dev_bufs[tid], roc_x_info[tid]);
      if(copy_from_dev)
      {
          hipMemcpyAsync(host_bufs[tid], dev_bufs[tid],  global_basis_alloc*sizeof(double), gpuMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          double *tptr1 = (double *)host_bufs[tid];
          if(out != tptr1) for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
#else
      if((double *)in == out)
          fftw_execute_dft_c2r (fftw_r2c_backward_plan_inplace,  reinterpret_cast<fftw_complex*>(in), out);
      else
          fftw_execute_dft_c2r (fftw_r2c_backward_plan,  reinterpret_cast<fftw_complex*>(in), out);
#endif
  }
  else
  {
  }
}

void Pw::FftInverse (std::complex<float> * in, float * out)
{
    FftInverse(in, out, true, true, true);
}

void Pw::FftInverse (std::complex<float> * in, float * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{

  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      if(copy_to_dev)
      {
          std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(std::complex<float>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecC2R(gpu_plans_c2r[tid], (cufftComplex*)dev_bufs[tid], (cufftReal*)dev_bufs[tid]);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
          gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid],  global_basis_alloc*sizeof(float), gpuMemcpyDeviceToHost, streams[tid]);
          gpuStreamSynchronize(streams[tid]);
          float *tptr1 = (float *)host_bufs[tid];
          if(out != tptr1) for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
#elif HIP_ENABLED
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
          for(size_t i = 0;i < global_basis_alloc;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], host_bufs[tid], global_basis_alloc*sizeof(std::complex<float>), gpuMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_c2r[tid], (void**) &dev_bufs[tid], (void**) &dev_bufs[tid], roc_x_info[tid]);
      if(copy_from_dev)
      {
          hipMemcpyAsync(host_bufs[tid], dev_bufs[tid],  global_basis_alloc*sizeof(float), gpuMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          float *tptr1 = (float *)host_bufs[tid];
          if(out != tptr1) for(size_t i = 0;i < global_basis_alloc;i++) out[i] = tptr1[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
#else
      if((float *)in == out)
          fftwf_execute_dft_c2r (fftwf_r2c_backward_plan_inplace,  reinterpret_cast<fftwf_complex*>(in), out);
      else
          fftwf_execute_dft_c2r (fftwf_r2c_backward_plan,  reinterpret_cast<fftwf_complex*>(in), out);
#endif
  }
  else
  {
  }
}

void Pw::FftInverse (std::complex<double> * in, std::complex<double> * out)
{
    FftInverse(in, out, true, true, true);
}

void Pw::FftInverse (std::complex<double> * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();

  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      if(copy_to_dev)
      {
	  for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
	  gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<double>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecZ2Z(gpu_plans[tid], (cufftDoubleComplex*)dev_bufs[tid], (cufftDoubleComplex*)dev_bufs[tid], CUFFT_INVERSE);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
	  gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<double>), gpuMemcpyDeviceToHost, streams[tid]);
	  gpuStreamSynchronize(streams[tid]);
	  for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
#elif HIP_ENABLED
      std::complex<double> *tptr = host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], tptr, pbasis*sizeof(std::complex<double>), hipMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_inv[tid], (void**) &dev_bufs[tid], NULL, roc_x_info[tid]);
      if(copy_from_dev)
      {   
          hipMemcpyAsync(tptr, dev_bufs[tid], pbasis*sizeof(std::complex<double>), hipMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
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
          for(size_t i = 0;i < pbasis;i++) buf[i] = in[i];
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
    FftInverse(in, out, true, true, true);
}

void Pw::FftInverse (std::complex<float> * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu)
{
  BaseThread *T = BaseThread::getBaseThread(0);
  int tid = T->get_thread_tid();
  if(tid < 0) tid = 0;
  if(tid == 0) tid = omp_get_thread_num();
  
  if(Grid->get_NPES() == 1)
  {
#if CUDA_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      if(copy_to_dev)
      {
	  for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
	  gpuMemcpyAsync(dev_bufs[tid], host_bufs[tid], pbasis*sizeof(std::complex<float>), gpuMemcpyHostToDevice, streams[tid]);
      }
      cufftExecC2C(gpu_plans_f[tid], (cufftComplex*)dev_bufs[tid], (cufftComplex*)dev_bufs[tid], CUFFT_INVERSE);
      gpuStreamSynchronize(streams[tid]);
      if(copy_from_dev)
      {
	  gpuMemcpyAsync(host_bufs[tid], dev_bufs[tid], pbasis*sizeof(std::complex<float>), gpuMemcpyDeviceToHost, streams[tid]);
	  gpuStreamSynchronize(streams[tid]);
	  for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
#elif HIP_ENABLED
      std::complex<float> *tptr = (std::complex<float> *)host_bufs[tid];
      hipStreamSynchronize(streams[tid]);
      if(copy_to_dev)
      {
          for(size_t i = 0;i < pbasis;i++) tptr[i] = in[i];
          hipMemcpyAsync(dev_bufs[tid], tptr, pbasis*sizeof(std::complex<float>), hipMemcpyHostToDevice, streams[tid]);
      }
      rocfft_execute(gpu_plans_f_inv[tid], (void**) &dev_bufs[tid], NULL, roc_x_info[tid]);
      if(copy_from_dev)
      {   
          hipMemcpyAsync(tptr, dev_bufs[tid], pbasis*sizeof(std::complex<float>), hipMemcpyDeviceToHost, streams[tid]);
          hipStreamSynchronize(streams[tid]);
          for(size_t i = 0;i < pbasis;i++) out[i] = tptr[i];
      }
      else
      {
          hipStreamSynchronize(streams[tid]);
      }
#else 
      if(in == out)
      {
          fftwf_execute_dft (fftwf_backward_plan_inplace,  reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
      }
      else
      {
          std::complex<float> *buf = new std::complex<float>[pbasis];
          for(size_t i = 0;i < pbasis;i++) buf[i] = in[i];
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
          for(size_t i = 0;i < pbasis;i++) buf[i] = in[i];
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

  int max_threads = std::max(ct.MG_THREADS_PER_NODE, ct.OMP_THREADS_PER_NODE);
  for(int thread=0;thread < max_threads;thread++)
  {
      if(this->distributed_plan[thread]) fft_3d_destroy_plan(this->distributed_plan[thread]);
      if(this->distributed_plan_f[thread]) fft_3d_destroy_plan(this->distributed_plan_f[thread]);
  }

  if(Grid->get_NPES() == 1)
  {
      fftw_destroy_plan(fftw_r2c_backward_plan_inplace);
      fftw_destroy_plan(fftw_r2c_forward_plan_inplace);
      fftw_destroy_plan(fftw_r2c_backward_plan);
      fftw_destroy_plan(fftw_r2c_forward_plan);

      fftw_destroy_plan(fftw_backward_plan_inplace);
      fftw_destroy_plan(fftw_forward_plan_inplace);
      fftw_destroy_plan(fftw_backward_plan);
      fftw_destroy_plan(fftw_forward_plan);

      fftwf_destroy_plan(fftwf_r2c_backward_plan_inplace);
      fftwf_destroy_plan(fftwf_r2c_forward_plan_inplace);
      fftwf_destroy_plan(fftwf_r2c_backward_plan);
      fftwf_destroy_plan(fftwf_r2c_forward_plan);

      fftwf_destroy_plan(fftwf_backward_plan_inplace);
      fftwf_destroy_plan(fftwf_forward_plan_inplace);
      fftwf_destroy_plan(fftwf_backward_plan);
      fftwf_destroy_plan(fftwf_forward_plan);
#if CUDA_ENABLED
      for (int i = 0; i < num_streams; i++)
      {
         cufftDestroy(gpu_plans[i]);
         cufftDestroy(gpu_plans_f[i]);
         cufftDestroy(gpu_plans_r2c[i]);
         cufftDestroy(gpu_plans_c2r[i]);
         cufftDestroy(gpu_plans_d2z[i]);
         cufftDestroy(gpu_plans_z2d[i]);
         RmgGpuError(__FILE__, __LINE__, gpuFreeHost(host_bufs[i]), "Error: gpuFreeHost failed.\n");
         RmgGpuError(__FILE__, __LINE__, gpuFree(dev_bufs[i]), "Error: gpuFreeHost failed.\n");
      }

      // Gpu streams and plans
      for (int i = 0; i < num_streams; i++)
          RmgGpuError(__FILE__, __LINE__, gpuStreamDestroy(streams[i]), "Problem freeing gpu stream.");


#endif

      int nthreads = std::max(ct.OMP_THREADS_PER_NODE, ct.MG_THREADS_PER_NODE);
      for (int i = 0; i < nthreads; i++)
      {
          fftw_free(host_bufs[i]);
      }
  }
}

