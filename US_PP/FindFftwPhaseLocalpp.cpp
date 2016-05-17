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
void FindFftwPhaseLocalpp (int nldim, double *nlcdrs, std::complex<double>* phase_fft, int grid_level)
{

    // grid_level = 1 for wavefunction grid
    // grid_level = 2 for rho_grid

    int i1, j1, k1;

    /*Reciprocal grid spacings in x, y and z directions */
    double rgs_x = 1.0 / (Rmg_G->get_hxgrid(grid_level) * Rmg_L.get_xside());
    double rgs_y = 1.0 / (Rmg_G->get_hygrid(grid_level) * Rmg_L.get_yside());
    double rgs_z = 1.0 / (Rmg_G->get_hzgrid(grid_level) * Rmg_L.get_zside());

    int nldim_sq = nldim * nldim;


    for (int i = -nldim / 2; i <= nldim / 2; i++)
    {
        for (int j = -nldim / 2; j <= nldim / 2; j++)
        {
            for (int k = -nldim / 2; k <= nldim / 2; k++)
            {

                if (i < 0)
                    i1 = i + nldim;
                else
                    i1 = i;

                if (j < 0)
                    j1 = j + nldim;
                else
                    j1 = j;

                if (k < 0)
                    k1 = k + nldim;
                else
                    k1 = k;


                /* Phase factor */
                double theta = 2.0 * PI / nldim *
                    (((nlcdrs[0] * (double) i) * rgs_x)
                     + ((nlcdrs[1] * (double) j) * rgs_y) + ((nlcdrs[2] * (double) k) * rgs_z));

                int idx1 = i1 * nldim_sq + j1 * nldim + k1;

                phase_fft[idx1] = exp(std::complex<double>(0.0, theta));
            }
        }
    }

}
