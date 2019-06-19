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


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"


/*This calculates the phase factor that will be used when calculating the backwards fourier transform*/
void FindPhaseKpoint (double *kvec, int nlxdim, int nlydim, int nlzdim, double * nlcdrs, std::complex<double>* phase_fftw, bool localize)
{

    int i1, j1, k1;
    double tpiba = 2.0*PI / Rmg_L.celldm[0];

    // If we are running with localized projectors then nxldim,nlydim and nlzdim
    // are odd numbers so the loop from (-nlxdim / 2) to (nlxdim / 2) is correct but
    // if running with non-localized projectors then they are odd and we need to adjust
    // the loop.
    int ixadj = 1;
    int iyadj = 1;
    int izadj = 1;
    if(nlxdim % 2) ixadj = 0;
    if(nlydim % 2) iyadj = 0;
    if(nlzdim % 2) izadj = 0;

    int ilo = 0;
    int jlo = 0;
    int klo = 0;
    int ihi = nlxdim;
    int jhi = nlydim;
    int khi = nlzdim;
    if(!localize) {
        ilo = get_PX_OFFSET();
        jlo = get_PY_OFFSET();
        klo = get_PZ_OFFSET();
        ihi = ilo + get_PX0_GRID();
        jhi = jlo + get_PY0_GRID();
        khi = klo + get_PZ0_GRID();
    } 

    for (int i = -nlxdim / 2; i <= nlxdim / 2 - ixadj; i++)
    {
        for (int j = -nlydim / 2; j <= nlydim / 2 - iyadj; j++)
        {
            for (int k = -nlzdim / 2; k <= nlzdim / 2 - izadj; k++)
            {

                if (i < 0)
                    i1 = i + nlxdim;
                else
                    i1 = i;

                if (j < 0)
                    j1 = j + nlydim;
                else
                    j1 = j;

                if (k < 0)
                    k1 = k + nlzdim;
                else
                    k1 = k;


                bool map = (i1 >= ilo) && (i1 < ihi) && (j1 >= jlo) && (j1 < jhi) && (k1 >= klo) && (k1 < khi);
                if(map) 
                {
                    /* Phase factor */
                    int idx1 = (i1 - ilo) * (jhi - jlo) * (khi - klo) + (j1 - jlo) * (khi - klo) + (k1 - klo);
                    double theta = nlcdrs[0] * (tpiba*coarse_pwaves->g[idx1].a[0] + kvec[0]) +
                                   nlcdrs[1] * (tpiba*coarse_pwaves->g[idx1].a[1] + kvec[1]) +
                                   nlcdrs[2] * (tpiba*coarse_pwaves->g[idx1].a[2] + kvec[2]);

                    phase_fftw[idx1] = exp(std::complex<double>(0.0, theta));
                }

            }
        }
    }

}
