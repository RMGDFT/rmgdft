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

void FftInitPlans(void);

void FftInterpolation (BaseGrid &G, double *coarse, double *fine, int ratio);

void FftGradientCoarse(double *x, double *fgx, double *fgy, double *fgz);

void FftGradientFine(double *x, double *fgx, double *fgy, double *fgz);

void FftGradient(double *x, double *fgx, double *fgy, double *fgz, Pw &pwaves);

void FftLaplacianCoarse(double *x, double *lapx);

void FftLaplacianFine(double *x, double *lapx);

void FftLaplacian(double *x, double *lapx, Pw &pwaves);

void FftFilter(double *x, Pw &pwaves, double factor, int type);

void FftFreqBin(double *x, Pw &pwaves, double *bins);

void FftRestrict(double *fine, double *coarse, int ratio);

void PfftForward(double *, std::complex<double> *, Pw &pwaves);

void PfftForward(std::complex<double> *, std::complex<double> *, Pw &pwaves);

void PfftInverse(std::complex<double> *, std::complex<double> *, Pw &pwaves);

void PfftInverse(std::complex<double> *, double *, Pw &pwaves);

#endif
#endif
