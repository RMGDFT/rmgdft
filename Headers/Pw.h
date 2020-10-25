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
*/


#ifndef RMG_Pw_H
#define RMG_Pw_H 1


// This class is used to handle the integration of plane wave functionality into RMG

#include <complex>
#include <vector>

#if CUDA_ENABLED
   #include <cuda.h>
   #include <cufft.h>
#endif

#include "BaseGrid.h"
#include "Lattice.h"
#include "fftw3.h"
#include "fft3d.h"
#include "WaitQueue.h"


typedef struct {
  double a[3];
} gvector;

typedef struct {
  std::complex<double> *in;
  std::complex<double> *out;
  int type;
} FftPair;

class Pw {

private:

    // Pending ffts
    WaitQueue<FftPair> FftQ;

    // Grid spacings
    double hxgrid;
    double hygrid;
    double hzgrid;

    // true if gamma only
    bool is_gamma;

    // Parallel fft plans
//    struct fft_plan_3d *distributed_plan;
    std::vector<struct fft_plan_3d<fftw_complex, double> *> distributed_plan;
    std::vector<struct fft_plan_3d<fftwf_complex, float> *> distributed_plan_f;

    // Local fft plans
    fftw_plan fftw_r2c_forward_plan, fftw_r2c_backward_plan;
    fftw_plan fftw_r2c_forward_plan_inplace, fftw_r2c_backward_plan_inplace;
    fftw_plan fftw_forward_plan, fftw_backward_plan;
    fftw_plan fftw_forward_plan_inplace, fftw_backward_plan_inplace;
    fftwf_plan fftwf_r2c_forward_plan, fftwf_r2c_backward_plan;
    fftwf_plan fftwf_r2c_forward_plan_inplace, fftwf_r2c_backward_plan_inplace;
    fftwf_plan fftwf_forward_plan, fftwf_backward_plan;
    fftwf_plan fftwf_forward_plan_inplace, fftwf_backward_plan_inplace;
 
    MPI_Comm comm;

public:
    Pw (BaseGrid &G, Lattice &L, int ratio, bool gamma_flag);
    void index_to_gvector(int *index, double *gvector);
    size_t count_filtered_gvectors(double filter_factor);
    void FftForward (double * in, std::complex<double> * out);
    void FftForward (std::complex<double> * in, std::complex<double> * out);
    void FftInverse (std::complex<double> * in, std::complex<double> * out);
    void FftInverse (std::complex<double> * in, double * out);

    void FftForward (float * in, std::complex<float> * out);
    void FftForward (std::complex<float> * in, std::complex<float> * out);
    void FftInverse (std::complex<float> * in, std::complex<float> * out);
    void FftInverse (std::complex<float> * in, float * out);

    void FftForward (double * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftForward (std::complex<double> * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftInverse (std::complex<double> * in, std::complex<double> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftInverse (std::complex<double> * in, double * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);

    void FftForward (float * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftForward (std::complex<float> * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftInverse (std::complex<float> * in, std::complex<float> * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);
    void FftInverse (std::complex<float> * in, float * out, bool copy_to_dev, bool copy_from_dev, bool use_gpu);

    ~Pw(void);

    // BaseGrid class
    BaseGrid *Grid;

    // Lattice object
    Lattice *L;

    // Real space basis on this node and globally
    size_t pbasis;
    size_t global_basis;

    // Packed global basis for EXX local ffts
    size_t global_basis_packed;

    // Real space grid dimensions on this node
    int dimx;
    int dimy;
    int dimz;

    // Global dimensions of the real space grid
    int global_dimx;
    int global_dimy;
    int global_dimz;

    // Plane wave cutoff and max
    double gcut;
    double gmax;

    size_t ng;
    gvector *g;
    double *gmags;
    bool *gmask;

    std::vector<std::complex<double> *> host_bufs;

#if CUDA_ENABLED
    int num_streams;
    std::vector<cudaStream_t> streams;
    std::vector<cufftHandle> gpu_plans;
    std::vector<cufftHandle> gpu_plans_f;
    std::vector<cufftHandle> gpu_plans_r2c;
    std::vector<cufftHandle> gpu_plans_c2r;
    std::vector<cufftHandle> gpu_plans_d2z;
    std::vector<cufftHandle> gpu_plans_z2d;
    std::vector<std::complex<double> *> dev_bufs;
#endif    

};


#endif
