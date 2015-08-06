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
#include <float.h>
#include <math.h>
#include "main.h"

void nlforce_par_omega (double * par_omega, int ion, int nh, double *force)
{
    int idx, idx1, size, n, m;
    double forces[3];
    double *omega_x, *omega_y, *omega_z, *qqq;

    size = nh * (nh + 1) / 2;

    omega_x = par_omega;
    omega_y = omega_x + size;
    omega_z = omega_y + size;

    qqq = pct.qqq[ion];

    for (idx = 0; idx < 3; idx++)
        forces[idx] = 0.0;

    idx = 0;
    for (n = 0; n < nh; n++)
    {
        for (m = n; m < nh; m++)
        {
            idx1 = n * nh + m;
            if (n == m)
            {
                forces[0] += qqq[idx1] * omega_x[idx];
                forces[1] += qqq[idx1] * omega_y[idx];
                forces[2] += qqq[idx1] * omega_z[idx];
            }
            else
            {
                forces[0] += 2.0 * qqq[idx1] * omega_x[idx];
                forces[1] += 2.0 * qqq[idx1] * omega_y[idx];
                forces[2] += 2.0 * qqq[idx1] * omega_z[idx];
            }

            ++idx;
        }
    }


    force[0] -= forces[0];
    force[1] -= forces[1];
    force[2] -= forces[2];

}
