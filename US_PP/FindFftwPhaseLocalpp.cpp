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
void FindFftwPhaseLocalpp (int nlxdim, int nlydim, int nlzdim, double *nlcdrs, std::complex<double>* phase_fft, int grid_level)
{

    // grid_level = 1 for wavefunction grid
    // grid_level = 2 for rho_grid

    int i1, j1, k1;

    /*Reciprocal grid spacings in x, y and z directions */
    double rgs_x = 2.0 * PI / (Rmg_G->get_hxgrid(grid_level) * Rmg_L.get_xside() * nlxdim);
    double rgs_y = 2.0 * PI / (Rmg_G->get_hygrid(grid_level) * Rmg_L.get_yside() * nlydim);
    double rgs_z = 2.0 * PI / (Rmg_G->get_hzgrid(grid_level) * Rmg_L.get_zside() * nlzdim);

    int nldim_sq = nlydim * nlzdim;


    for (int i = -nlxdim / 2; i <= nlxdim / 2; i++)
    {
        for (int j = -nlzdim / 2; j <= nlydim / 2; j++)
        {
            for (int k = -nlzdim / 2; k <= nlzdim / 2; k++)
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


                /* Phase factor */
                double theta = ((nlcdrs[0] * (double) i) * rgs_x)
                     + ((nlcdrs[1] * (double) j) * rgs_y) + ((nlcdrs[2] * (double) k) * rgs_z);

                int idx1 = i1 * nldim_sq + j1 * nlzdim + k1;

                phase_fft[idx1] = exp(std::complex<double>(0.0, theta));
            }
        }
    }

}
