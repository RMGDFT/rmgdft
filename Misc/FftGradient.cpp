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
#include "vdW.h"
#include "RmgException.h"
#include "xc.h"
#include "RmgSumAll.h"
#include "transition.h"

#if USE_PFFT
#include "RmgParralelFft.h"

void FftGradient(double *x, double *fgx, double *fgy, double *fgz, Pw &pwaves)
{

    pfft_plan forw, inv;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = 1.0 / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    ptrdiff_t grid[3]; 
    grid[0] = pwaves.global_dimx;
    grid[1] = pwaves.global_dimy;
    grid[2] = pwaves.global_dimz;

    std::complex<double> czero(0.0,0.0);
    std::complex<double> ci(0.0,1.0);
    std::complex<double> *tx = new std::complex<double>[isize];
    std::complex<double> *cgx = new std::complex<double>[isize];

    forw = pfft_plan_dft_3d(grid,
                          (double (*)[2])tx,
                          (double (*)[2])tx,
                          pct.pfft_comm,
                          PFFT_FORWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    inv = pfft_plan_dft_3d(grid,
                          (double (*)[2])cgx,
                          (double (*)[2])cgx,
                          pct.pfft_comm,
                          PFFT_BACKWARD,
                          PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = std::complex<double>(x[ix], 0.0);
    }

    pfft_execute_dft(forw, (double (*)[2])tx, (double (*)[2])tx);

    for(int icar=0;icar < 3;icar++) {

        for(int ix=0;ix < isize;ix++) cgx[ix] = czero;

        for(int ig=0;ig < isize;ig++) {
            if(pwaves.gmask[ig] == 1.0) {
                cgx[ig] = ci * tpiba * pwaves.g[ig].a[icar] * tx[ig];
            }
            else {
                cgx[ig] = czero;
            }
        }

        if(pct.gridpe == 0) cgx[0] = czero;
        pfft_execute_dft(inv, (double (*)[2])cgx, (double (*)[2])cgx);

        double *ts;
        if(icar == 0) ts = fgx;
        if(icar == 1) ts = fgy;
        if(icar == 2) ts = fgz;

        for(int ix=0;ix < isize;ix++) ts[ix] = scale * std::real(cgx[ix]);
    }

    pfft_destroy_plan(forw);
    pfft_destroy_plan(inv);

    delete [] cgx;
    delete [] tx;
}

#endif
