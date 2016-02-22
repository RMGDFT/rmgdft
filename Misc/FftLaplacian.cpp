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

#if USE_PFFT
#include "RmgParallelFft.h"

void FftLaplacianCoarse(double *x, double *lapx)
{
    FftLaplacian(x, lapx, *coarse_pwaves);
}

void FftLaplacianFine(double *x, double *lapx)
{
    FftLaplacian(x, lapx, *fine_pwaves);
}

void FftLaplacian(double *x, double *lapx, Pw &pwaves)
{

    pfft_plan forw, inv;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<double> czero(0.0,0.0);
    std::complex<double> *tx = new std::complex<double>[isize];

    if(&pwaves == coarse_pwaves) {
        forw = forward_coarse;
        inv = backward_coarse;
    }
    else if(&pwaves == fine_pwaves) {
        forw = forward_fine;
        inv = backward_fine;
    }
    else {
        ptrdiff_t grid[3]; 
        grid[0] = pwaves.global_dimx;
        grid[1] = pwaves.global_dimy;
        grid[2] = pwaves.global_dimz;

        forw = pfft_plan_dft_3d(grid,
                              (double (*)[2])tx,
                              (double (*)[2])tx,
                              pct.pfft_comm,
                              PFFT_FORWARD,
                              PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

        inv = pfft_plan_dft_3d(grid,
                              (double (*)[2])tx,
                              (double (*)[2])tx,
                              pct.pfft_comm,
                              PFFT_BACKWARD,
                              PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);
    }

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = std::complex<double>(x[ix], 0.0);
    }

    pfft_execute_dft(forw, (double (*)[2])tx, (double (*)[2])tx);

    for(int ig=0;ig < isize;ig++) {
        if(pwaves.gmask[ig] == 1.0) {
            tx[ig] = -pwaves.gmags[ig] * tx[ig];
        }
        else {
            tx[ig] = czero;
        }
    }
    if(pct.gridpe == 0) tx[0] = czero;

    pfft_execute_dft(inv, (double (*)[2])tx, (double (*)[2])tx);

    for(int ix=0;ix < isize;ix++) lapx[ix] = scale * std::real(tx[ix]);

    if((&pwaves != coarse_pwaves) && (&pwaves != fine_pwaves)) {
        pfft_destroy_plan(inv);
        pfft_destroy_plan(forw);
    }

    delete [] tx;
}

#endif
