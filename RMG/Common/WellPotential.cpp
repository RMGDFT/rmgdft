/*
 *
 * Copyright 2024 The RMG Project Developers. See the COPYRIGHT file 
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

#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"

// Experimental code to implement a well potential for non-periodic boundary conditions.

void WellPotential (double *v)
{
    double vector[3];
    double x[3], xtal[3];
    xtal[0] = 0.5;
    xtal[1] = 0.5;
    xtal[2] = 0.5;

    int ratio = Rmg_G->default_FG_RATIO;
    double hxgrid = Rmg_G->get_hxgrid(ratio);
    double hygrid = Rmg_G->get_hygrid(ratio);
    double hzgrid = Rmg_G->get_hzgrid(ratio);
    int dimx = Rmg_G->get_PX0_GRID(ratio);
    int dimy = Rmg_G->get_PY0_GRID(ratio);
    int dimz = Rmg_G->get_PZ0_GRID(ratio);
    int ilow = Rmg_G->get_PX_OFFSET(ratio);
    int jlow = Rmg_G->get_PY_OFFSET(ratio);
    int klow = Rmg_G->get_PZ_OFFSET(ratio);

    int ihi = ilow + dimx;
    int jhi = jlow + dimy;
    int khi = klow + dimz;

    // r will hold the distance from the center of the cell to a specific grid point
    std::vector<double> r;
    r.resize(dimx*dimy*dimz);
    double rmax = 0.0;
    for (int ix = ilow; ix < ihi; ix++)
    {
        x[0] = ix * hxgrid - xtal[0];
        for (int iy = jlow; iy < jhi; iy++)
        {
            x[1] = iy * hygrid - xtal[1];
            for (int iz = klow; iz < khi; iz++)
            {
                x[2] = iz * hzgrid - xtal[2];

                int idx = (ix-ilow)*dimy*dimz + (iy-jlow)*dimz + iz-klow;
                r[idx] = Rmg_L.metric(x);
                rmax = std::max(r[idx], rmax);
//                Rmg_L.to_cartesian(x, vector);

            }
        }
    }
    
    // Get max over all PEs
    MPI_Allreduce(MPI_IN_PLACE, &rmax, 1, MPI_DOUBLE, MPI_MAX, pct.grid_comm);

    // Confinement properties of the well depend on it's height and width.
    // Units of vmax are hartrees so 10 corresponds to a height of 272 eV.
    double vmax = 10.0;

    // Fermi distribution
    rmax = rmax - 5.0*hxgrid*Rmg_L.get_xside();
    double w = hxgrid*Rmg_L.get_xside()/5.0;
    for (int idx = 0; idx < dimx*dimy*dimz;idx++)
    {
        double fd = 1.0 - 1.0 / (1.0 + exp((r[idx] - rmax)/w));
        double v1 = vmax*fd;
        v[idx] += v1;
    }

}
