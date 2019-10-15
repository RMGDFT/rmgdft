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


// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)
#ifdef FFT_SINGLE
#define FFT_PRECISION 1
#else
#define FFT_PRECISION 2
#endif

#define FFT_FFTW 1
// set default fftw library. switch to FFT_FFTW3 when convenient.
#ifdef FFT_FFTW
#define FFT_FFTW3
#endif

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION == 1

#define FFTW_API(function)  fftwf_ ## function

// -------------------------------------------------------------------------

// Data types for double-precision complex

#elif FFT_PRECISION == 2

#define FFTW_API(function)  fftw_ ## function

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
  FFTW_API(plan) plan_fast_forward;
  FFTW_API(plan) plan_fast_backward;
  FFTW_API(plan) plan_mid_forward;
  FFTW_API(plan) plan_mid_backward;
  FFTW_API(plan) plan_slow_forward;
  FFTW_API(plan) plan_slow_backward;
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

/* ERROR/WARNING messages:

*/
#endif
