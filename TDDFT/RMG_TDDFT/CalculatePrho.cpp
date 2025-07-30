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



void CalculatePrho (double *rho, double *prho_t)
{

    int pbasis = fine_pwaves->pbasis;
    int size = pbasis;

    std::complex<double> *a_in = new std::complex<double>[size];
    std::complex<double> *a_out = new std::complex<double>[size];


    for(int i = 0;i < pbasis;i++) a_in[i] = rho[i];

    fine_pwaves->FftForward(a_in, a_out);

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    for(int ig=0;ig < pbasis;ig++) {
        if((fine_pwaves->gmags[ig] > 1.0e-6) && fine_pwaves->gmask[ig]) 
            a_in[ig] = a_out[ig]/(fine_pwaves->gmags[ig] *tpiba2) * (fine_pwaves->g[ig].a[0] * tpiba);
        else
            a_in[ig] = 0.0;
    }

    fine_pwaves->FftInverse(a_in, a_in);
    prho_t[0] = 0.0;
    for(int i = 0;i < pbasis;i++) prho_t[0] += std::real(a_in[i])/(double)fine_pwaves->global_basis;
    
    for(int ig=0;ig < pbasis;ig++) {
        if((fine_pwaves->gmags[ig] > 1.0e-6) && fine_pwaves->gmask[ig]) 
            a_in[ig] = a_out[ig]/(fine_pwaves->gmags[ig] *tpiba2) * (fine_pwaves->g[ig].a[1] * tpiba);
        else
            a_in[ig] = 0.0;
    }

    fine_pwaves->FftInverse(a_in, a_in);
    prho_t[1] = 0.0;
    for(int i = 0;i < pbasis;i++) prho_t[1] += std::real(a_in[i])/(double)fine_pwaves->global_basis;
    
    for(int ig=0;ig < pbasis;ig++) {
        if((fine_pwaves->gmags[ig] > 1.0e-6) && fine_pwaves->gmask[ig]) 
            a_in[ig] = a_out[ig]/(fine_pwaves->gmags[ig] *tpiba2) * (fine_pwaves->g[ig].a[2] * tpiba);
        else
            a_in[ig] = 0.0;
    }

    fine_pwaves->FftInverse(a_in, a_in);
    prho_t[2] = 0.0;
    for(int i = 0;i < pbasis;i++) prho_t[2] += std::real(a_in[i])/(double)fine_pwaves->global_basis;

    MPI_Allreduce(MPI_IN_PLACE, prho_t, 3, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    delete [] a_in;
    delete [] a_out;

}

