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

void VhPfft(double *rho_tot, double *rhoc, double *vh)
{

    ptrdiff_t densgrid[3];
    pfft_plan plan_forward, plan_back;
    densgrid[0] = fine_pwaves->global_dimx;
    densgrid[1] = fine_pwaves->global_dimy;
    densgrid[2] = fine_pwaves->global_dimz;

    int global_basis = fine_pwaves->global_basis;
    int pbasis = fine_pwaves->pbasis;

    std::complex<double> *crho = new std::complex<double>[pbasis];
    plan_forward = pfft_plan_dft_3d(densgrid,
            (double (*)[2])crho,
            (double (*)[2])crho,
            pct.pfft_comm,
            PFFT_FORWARD,
            PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    plan_back = pfft_plan_dft_3d(densgrid,
            (pfft_complex *)crho,
            (pfft_complex *)crho,
            pct.pfft_comm,
            PFFT_BACKWARD,
            PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    for(int i = 0;i < pbasis;i++) crho[i] = std::complex<double>(rho_tot[i]-rhoc[i], 0.0);
    pfft_execute_dft(plan_forward, (double (*)[2])crho, (double (*)[2])crho);

    double tem = 0.0;
    for(int i = 0;i < pbasis;i++) tem += rho_tot[i] - rhoc[i];
    if(pct.gridpe == 0 && abs(crho[0]) >1.0e-6) 
        printf("\n WARNING:  total charge is not zero: crho = %e %e rho-rhoc = %e", std::real(crho[0]), std::imag(crho[0]), tem); 

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) crho[ig] = crho[ig]/(fine_pwaves->gmags[ig] *tpiba2);
    }

    pfft_execute_dft(plan_back, (double (*)[2])crho, (double (*)[2])crho);
    for(int i = 0;i < pbasis;i++) vh[i] = std::real(crho[i])/(double)global_basis * 4.0 * PI;

    pfft_destroy_plan(plan_back);
    pfft_destroy_plan(plan_forward);
    delete [] crho;

}

