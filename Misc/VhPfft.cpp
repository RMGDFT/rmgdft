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

void VhPfft(double *rho_tot, double *rhoc, double *vh)
{

    int pbasis = fine_pwaves->pbasis;
    int size = pbasis;
    std::complex<double> ZERO_t(0.0, 0.0);
    std::complex<double> *crho = new std::complex<double>[size];


    for(int i = 0;i < pbasis;i++) crho[i] = std::complex<double>(rho_tot[i], 0.0);
    fine_pwaves->FftForward(crho, crho);

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    for(int ig=0;ig < pbasis;ig++) {
        if((fine_pwaves->gmags[ig] > 1.0e-6) && fine_pwaves->gmask[ig]) 
            crho[ig] = crho[ig]/(fine_pwaves->gmags[ig] *tpiba2);
        else
            crho[ig] = ZERO_t;
    }

    fine_pwaves->FftInverse(crho, crho);
    for(int i = 0;i < pbasis;i++) vh[i] = std::real(crho[i])/(double)fine_pwaves->global_basis * 4.0 * PI;

    delete [] crho;

}

