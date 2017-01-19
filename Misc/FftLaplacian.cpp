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

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double scale = tpiba * tpiba / (double)pwaves.global_basis;
    int isize = pwaves.pbasis;

    std::complex<double> czero(0.0,0.0);
    std::complex<double> *tx = new std::complex<double>[isize];

    for(int ix = 0;ix < isize;ix++) {
        tx[ix] = std::complex<double>(x[ix], 0.0);
    }

    PfftForward(tx, tx, pwaves);

    for(int ig=0;ig < isize;ig++) {
        if(pwaves.gmask[ig] == 1.0) {
            tx[ig] = -pwaves.gmags[ig] * tx[ig];
        }
        else {
            tx[ig] = czero;
        }
    }
    if(pct.gridpe == 0) tx[0] = czero;

    PfftInverse(tx, tx, pwaves);

    for(int ix=0;ix < isize;ix++) lapx[ix] = scale * std::real(tx[ix]);


    delete [] tx;
}

