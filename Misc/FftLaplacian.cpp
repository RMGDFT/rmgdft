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
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "RmgParallelFft.h"



void FftLaplacianCoarse(float *x, float *lapx)
{
   FftLaplacian(x, lapx, *coarse_pwaves);
}

void FftLaplacianCoarse(std::complex<float> *x, std::complex<float> *lapx)
{
   FftLaplacian(x, lapx, *coarse_pwaves);
}

void FftLaplacianCoarse(double *x, double *lapx)
{
    FftLaplacian(x, lapx, *coarse_pwaves);
}

void FftLaplacianCoarse(std::complex<double> *x, std::complex<double> *lapx)
{
    FftLaplacian(x, lapx, *coarse_pwaves);
}

void FftLaplacianFine(float *x, float *lapx)
{
   throw RmgFatalException() << "Float version not implemented yet in FftLaplacian "<< " at line " << __LINE__ << "\n";
}

void FftLaplacianFine(std::complex<float> *x, std::complex<float> *lapx)
{
   throw RmgFatalException() << "Float version not implemented yet in FftLaplacian "<< " at line " << __LINE__ << "\n";
}

void FftLaplacianFine(double *x, double *lapx)
{
    FftLaplacian(x, lapx, *fine_pwaves);
}

void FftLaplacianFine(std::complex<double> *x, std::complex<double> *lapx)
{
    FftLaplacian(x, lapx, *fine_pwaves);
}

void FftLaplacian(float *x, float *lapx, Pw &pwaves)
{

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<float> *tx = new std::complex<float>[isize];

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = std::complex<float>(x[ix], 0.0);
    }

    pwaves.FftForward(tx, tx);

    for(int ig=0;ig < isize;ig++) tx[ig] = -tx[ig] * (float)pwaves.gmags[ig];

    pwaves.FftInverse(tx, tx);

    for(int ix=0;ix < isize;ix++) lapx[ix] = (float)scale * (float)std::real(tx[ix]);
    delete [] tx;

}

void FftLaplacian(double *x, double *lapx, Pw &pwaves)
{

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<double> czero(0.0,0.0);
    std::complex<double> *tx = new std::complex<double>[isize];

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = std::complex<double>(x[ix], 0.0);
    }

    pwaves.FftForward(tx, tx);

    for(int ig=0;ig < isize;ig++) tx[ig] = -pwaves.gmags[ig] * tx[ig];

    pwaves.FftInverse(tx, tx);

    for(int ix=0;ix < isize;ix++) lapx[ix] = scale * std::real(tx[ix]);


    delete [] tx;
}

void FftLaplacian(std::complex<float> *x, std::complex<float> *lapx, Pw &pwaves)
{

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<float> *tx = new std::complex<float>[isize];

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = x[ix];
    }

    pwaves.FftForward(tx, tx);

    for(int ig=0;ig < isize;ig++) tx[ig] = -(float)pwaves.gmags[ig] * tx[ig];

    pwaves.FftInverse(tx, tx);

    for(int ix=0;ix < isize;ix++) lapx[ix] = (float)scale * tx[ix];


    delete [] tx;
}

void FftLaplacian(std::complex<double> *x, std::complex<double> *lapx, Pw &pwaves)
{

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<double> *tx = new std::complex<double>[isize];

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = x[ix];
    }

    pwaves.FftForward(tx, tx);

    for(int ig=0;ig < isize;ig++) tx[ig] = -pwaves.gmags[ig] * tx[ig];

    pwaves.FftInverse(tx, tx);

    for(int ix=0;ix < isize;ix++) lapx[ix] = scale * tx[ix];


    delete [] tx;
}

