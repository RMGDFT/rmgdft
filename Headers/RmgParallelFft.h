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

#ifndef RMG_RmgParallelFft_H
#define RMG_RmgParallelFft_H 1

#if __cplusplus
#include <complex>
#include "BaseGrid.h"
#include "Lattice.h"
#include "Pw.h"

typedef struct
{
    fftw_plan cplan;
    int nx;
    int ny;
    int nz;
    int pbasis;
#if RMG_CUDA_ENABLED
    cufftHandle gplan;
    cufftType type;
    int batch;
#endif
    
} LocalFftPlan;

void FftInitPlans(void);

void FftInterpolation (BaseGrid &G, double *coarse, double *fine, int ratio, bool use_sqrt);

void FftGradientCoarse(float *x, float *fgx, float *fgy, float *fgz);
void FftGradientCoarse(std::complex<float> *x, std::complex<float> *fgx, std::complex<float> *fgy, std::complex<float> *fgz);
void FftGradientCoarse(double *x, double *fgx, double *fgy, double *fgz);
void FftGradientCoarse(std::complex<double> *x, std::complex<double> *fgx, std::complex<double> *fgy, std::complex<double> *fgz);

void FftGradientFine(float *x, float *fgx, float *fgy, float *fgz);
void FftGradientFine(std::complex<float> *x, std::complex<float> *fgx, std::complex<float> *fgy, std::complex<float> *fgz);
void FftGradientFine(double *x, double *fgx, double *fgy, double *fgz);
void FftGradientFine(std::complex<double> *x, std::complex<double> *fgx, std::complex<double> *fgy, std::complex<double> *fgz);

void FftGradient(float *x, float *fgx, float *fgy, float *fgz, Pw &pwaves);
void FftGradient(double *x, double *fgx, double *fgy, double *fgz, Pw &pwaves);
void FftGradient(std::complex<float> *x, std::complex<float> *fgx, std::complex<float> *fgy, std::complex<float> *fgz, Pw &pwaves);
void FftGradient(std::complex<double> *x, std::complex<double> *fgx, std::complex<double> *fgy, std::complex<double> *fgz, Pw &pwaves);


void FftLaplacianCoarse(float *x, float *lapx);
void FftLaplacianCoarse(std::complex<float> *x, std::complex<float> *lapx);
void FftLaplacianCoarse(double *x, double *lapx);
void FftLaplacianCoarse(std::complex<double> *x, std::complex<double> *lapx);

void FftLaplacianFine(float *x, float *lapx);
void FftLaplacianFine(std::complex<float> *x, std::complex<float> *lapx);
void FftLaplacianFine(double *x, double *lapx);
void FftLaplacianFine(std::complex<double> *x, std::complex<double> *lapx);

void FftLaplacian(float *x, float *lapx, Pw &pwaves);
void FftLaplacian(std::complex<float> *x, std::complex<float> *lapx, Pw &pwaves);
void FftLaplacian(double *x, double *lapx, Pw &pwaves);
void FftLaplacian(std::complex<double> *x, std::complex<double> *lapx, Pw &pwaves);

void FftFilter(double *x, Pw &pwaves, Pw &c_pwaves, int type);
void FftFilter(double *x, Pw &pwaves, Pw &c_pwaves, double factor);

void FftFreqBin(double *x, Pw &pwaves, double *bins);

void FftRestrict(double *fine, double *coarse, int ratio);

void LocalFftForward(double *, std::complex<double> *, Pw &pwaves);

void LocalFftForward(std::complex<double> *, std::complex<double> *, Pw &pwaves);

void LocalFftInverse(std::complex<double> *, std::complex<double> *, Pw &pwaves);

void LocalFftInverse(std::complex<double> *, double *, Pw &pwaves);

void FftSmoother(double *x, Pw &pwaves, double factor);

#endif
#endif
