/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef RMG_fft3d_H
#define RMG_fft3d_H 1

#include <type_traits>
#include <complex>
#include "remap.h"
#include "fftw3.h"


#define FFT_FFTW 1
// set default fftw library. switch to FFT_FFTW3 when convenient.
#ifdef FFT_FFTW
#define FFT_FFTW3
#endif



// -------------------------------------------------------------------------

// details of how to do a 3d FFT

template <typename FFT_DATA, typename FFT_SCALAR>
struct fft_plan_3d {
  struct remap_plan_3d<FFT_SCALAR> *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_3d<FFT_SCALAR> *mid1_plan;      // remap from 1st -> 2nd FFTs
  struct remap_plan_3d<FFT_SCALAR> *mid2_plan;      // remap from 2nd -> 3rd FFTs
  struct remap_plan_3d<FFT_SCALAR> *post_plan;      // remap from 3rd FFTs -> output
  FFT_DATA *copy;                   // memory for remap results (if needed)
  FFT_DATA *scratch;                // scratch space for remaps
  int total1,total2,total3;         // # of 1st,2nd,3rd FFTs (times length)
  int length1,length2,length3;      // length of 1st,2nd,3rd FFTs
  int pre_target;                   // where to put remap results
  int mid1_target,mid2_target;
  int scaled;                       // whether to scale FFT results
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

                                    // system specific 1d FFT info
  fftw_plan plan_fast_forward;
  fftw_plan plan_fast_backward;
  fftw_plan plan_mid_forward;
  fftw_plan plan_mid_backward;
  fftw_plan plan_slow_forward;
  fftw_plan plan_slow_backward;
};

template<> struct fft_plan_3d<fftwf_complex,float>;

template <>
struct fft_plan_3d<fftwf_complex,float> {
  struct remap_plan_3d<float> *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_3d<float> *mid1_plan;      // remap from 1st -> 2nd FFTs
  struct remap_plan_3d<float> *mid2_plan;      // remap from 2nd -> 3rd FFTs
  struct remap_plan_3d<float> *post_plan;      // remap from 3rd FFTs -> output
  fftwf_complex *copy;                   // memory for remap results (if needed)
  fftwf_complex *scratch;                // scratch space for remaps
  int total1,total2,total3;         // # of 1st,2nd,3rd FFTs (times length)
  int length1,length2,length3;      // length of 1st,2nd,3rd FFTs
  int pre_target;                   // where to put remap results
  int mid1_target,mid2_target;
  int scaled;                       // whether to scale FFT results
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

                                    // system specific 1d FFT info
  fftwf_plan plan_fast_forward;
  fftwf_plan plan_fast_backward;
  fftwf_plan plan_mid_forward;
  fftwf_plan plan_mid_backward;
  fftwf_plan plan_slow_forward;
  fftwf_plan plan_slow_backward;
};



// function prototypes

template<typename FFT_DATA, typename FFT_SCALAR> void fft_3d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d<FFT_DATA, FFT_SCALAR> *);
template<typename FFT_DATA, typename FFT_SCALAR> struct fft_plan_3d<FFT_DATA, FFT_SCALAR> *fft_3d_create_plan(MPI_Comm, int, int, int,
                                       int, int, int, int, int,
                                       int, int, int, int, int, int, int,
                                       int, int, int *, int);
template<typename FFT_DATA, typename FFT_SCALAR> void fft_3d_destroy_plan(struct fft_plan_3d<FFT_DATA, FFT_SCALAR> *);
void factor(int, int *, int *);
void bifactor(int, int *, int *);
template<typename FFT_DATA, typename FFT_SCALAR> void fft_1d_only(FFT_DATA *, int, int, struct fft_plan_3d<FFT_DATA, FFT_SCALAR> *);

void fft_destroy_plan_wrapper(fftw_plan plan);
void fft_destroy_plan_wrapper(fftwf_plan plan);

void fftw_execute_plan_wrapper(fftw_plan plan, fftw_complex *d1, fftw_complex *d2);
void fftw_execute_plan_wrapper(fftwf_plan plan, fftwf_complex *d1, fftwf_complex *d2);

fftw_plan fft_plan_many_dft_wrapper(int rank, const int *n, int howmany,
                             fftw_complex *in, const int *inembed,
                             int istride, int idist,
                             fftw_complex *out, const int *onembed,
                             int ostride, int odist,
                             int sign, unsigned flags);

fftwf_plan fft_plan_many_dft_wrapper(int rank, const int *n, int howmany,
                             fftwf_complex *in, const int *inembed,
                             int istride, int idist,
                             fftwf_complex *out, const int *onembed,
                             int ostride, int odist,
                             int sign, unsigned flags);

#endif
