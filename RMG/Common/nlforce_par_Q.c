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
    int idx2, n, m, count, icount, size;
    int *pidx;

    count = pct.Qidxptrlen[ion];
    pidx = pct.Qindex[ion];
    double *Qnm;

    Qnm = pct.augfunc[ion];

    if (count)
    {
        size = (nh * (nh + 1)) / 2;

        idx2 = 0;
        for (n = 0; n < nh; n++)
        {
            for (m = n; m < nh; m++)
            {
                if (m != n) gamma[idx2] *= 2.0;
                idx2++;
            }
        }

        double one = 1.0, zero = 0.0, *tmp_arr;
        int ione = 1;
        tmp_arr = (double *)malloc(count * sizeof(double));

        dgemm("N", "N", &count, &ione, &size, &one, Qnm, &count, gamma, &size, &zero, tmp_arr, &count); 

        for (icount = 0; icount < count; icount++)
        {
            forces[0] -= gx[pidx[icount]] * tmp_arr[icount];
            forces[1] -= gy[pidx[icount]] * tmp_arr[icount];
            forces[2] -= gz[pidx[icount]] * tmp_arr[icount];
        }
        free(tmp_arr);
    }


}
