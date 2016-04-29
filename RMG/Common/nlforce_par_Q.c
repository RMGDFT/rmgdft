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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"

void nlforce_par_Q (double *gx, double *gy, double *gz, double * gamma, int ion, ION * iptr, int nh, double * forces)
{
    int idx1, idx2, n, m, count, icount, size;
    int *pidx;
    double tmp;


    /*Forces array is assumed to be already initialized */
    /*for(idx1=0;idx1<3;idx1++)  forces[idx1]=0.0; */


    count = pct.Qidxptrlen[ion];
    pidx = pct.Qindex[ion];
    double *Qnm;

    Qnm = pct.augfunc[ion];

    if (count)
    {
        size = (nh * (nh + 1) / 2) * count;


        for (icount = 0; icount < count; icount++)
        {
            tmp = 0.0;

            idx2 = 0;
            for (n = 0; n < nh; n++)
            {
                for (m = n; m < nh; m++)
                {
                    idx1 = idx2 * count + icount;
                    if (m == n)
                    {
                        tmp += Qnm[idx1] * gamma[idx2];
                    }
                    else
                    {
                        tmp += 2.0 * Qnm[idx1] * gamma[idx2];
                    }

                    ++idx2;
                }
            }
            forces[0] -= gx[pidx[icount]] * tmp;
            forces[1] -= gy[pidx[icount]] * tmp;
            forces[2] -= gz[pidx[icount]] * tmp;
        }

    }


}
