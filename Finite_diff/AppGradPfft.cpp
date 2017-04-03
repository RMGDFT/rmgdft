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
#include "RmgSumAll.h"
#include "transition.h"
#include "RmgParallelFft.h"
#include "RmgException.h"
#include "Pw.h"



template void AppGradPfft(double *a, double *gx, double *gy, double *gz, char *grid);
template void AppGradPfft(float *a, float *gx, float *gy, float *gz, char *grid);
template void AppGradPfft(std::complex<float> *a, std::complex<float> *gx, std::complex<float> *gy, std::complex<float> *gz, char *grid);
template void AppGradPfft(std::complex<double> *a, std::complex<double> *gx, std::complex<double> *gy, std::complex<double> *gz, char *grid);
template <typename DataType>
void AppGradPfft (DataType *a, DataType *gx, DataType *gy, DataType *gz, char *grid)
{

    Pw *pwaves;

    const char *coarse = "Coarse";
    const char *fine = "Fine";

    if(!strcmp(grid, coarse)) {
        pwaves = coarse_pwaves;
    }
    else if(!strcmp(grid, fine)) {
        pwaves = fine_pwaves;
    }
    else {
        throw RmgFatalException() << "Error! Grid type " << grid << " not defined in "
                                 << __FILE__ << " at line " << __LINE__ << "\n";
    }

    
    int pbasis = pwaves->pbasis;
    int size = pbasis;

    std::complex<double> *a_in = new std::complex<double>[size];
    std::complex<double> *a_out = new std::complex<double>[size];


    for(int i = 0;i < pbasis;i++) a_in[i] = a[i];

    PfftForward(a_in, a_out, *pwaves);

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];

//    for(int ig=0;ig < pbasis;ig++) {
//        if(pwaves->gmags[ig] > 0.8 * pwaves->gcut) a_out[ig] = 0.0;
//   }

    for(int ig=0;ig < pbasis;ig++) {
        a_in[ig] = a_out[ig] * std::complex<double>(0.0, pwaves->g[ig].a[0] * tpiba);

    }
    PfftInverse(a_in, a_in, *pwaves);
    for(int i = 0;i < pbasis;i++) gx[i] = std::real(a_in[i])/(double)pwaves->global_basis;

    for(int ig=0;ig < pbasis;ig++) {
        a_in[ig] = a_out[ig] * std::complex<double>(0.0, pwaves->g[ig].a[1] * tpiba);

    }
    PfftInverse(a_in, a_in, *pwaves);
    for(int i = 0;i < pbasis;i++) gy[i] = std::real(a_in[i])/(double)pwaves->global_basis;

    for(int ig=0;ig < pbasis;ig++) {
        a_in[ig] = a_out[ig] * std::complex<double>(0.0, pwaves->g[ig].a[2] * tpiba);

    }
    PfftInverse(a_in, a_in, *pwaves);
    for(int i = 0;i < pbasis;i++) gz[i] = std::real(a_in[i])/(double)pwaves->global_basis;

    delete [] a_in;
    delete [] a_out;

}

