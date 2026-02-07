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
#include "main.h"
#include "transition.h"
#include "prototypes_on.h"

void write_rho_x(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;
    std::vector<double> zvec(get_FNX_GRID(), 0.0);
    /* Get this processors offset */
    poff = get_FPX_OFFSET();


    /* Loop over this processor */
    for (ix = 0; ix < get_FPX0_GRID(); ix++)
    {
        t1 = 0.0;
        for (iy = 0; iy < get_FPY0_GRID(); iy++)
            for (iz = 0; iz < get_FPZ0_GRID(); iz++)

                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[ix + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNX_GRID();
    rmg::reduce(zvec.data(), ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        rmg::printlog("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < get_FNX_GRID(); ix++)
        {
            t1 = ix * get_hxxgrid();
            rmg::printlog(" %d %f %s\n", ix, zvec[ix] / get_FNY_GRID() / get_FNZ_GRID(), ab);
        }
        rmg::printlog(" & %s\n", ab);
        fflush(NULL);
    }

}                               /* end get_avgd */

void write_rho_y(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;

    std::vector<double> zvec(get_FNY_GRID(), 0.0);
    /* Get this processors offset */
    poff = get_FPY_OFFSET();

    /* Loop over this processor */
    for (iy = 0; iy < get_FPY0_GRID(); iy++)
    {
        t1 = 0.0;
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
            for (iz = 0; iz < get_FPZ0_GRID(); iz++)

                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[iy + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNY_GRID();
    rmg::reduce(zvec.data(), ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        rmg::printlog("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < get_FNY_GRID(); ix++)
        {
            t1 = ix * get_hxxgrid();
            rmg::printlog(" %d %f %s\n", ix, zvec[ix] / get_FNX_GRID() / get_FNZ_GRID(), ab);
        }
        rmg::printlog(" & %s\n", ab);
        fflush(NULL);
    }

}                               /* end get_avgd */

void write_rho_z(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;

    std::vector<double> zvec(get_FNZ_GRID(), 0.0);
    /* Get this processors offset */
    poff = get_FPZ_OFFSET();
    int pxoff = get_FPX_OFFSET();
    int pyoff = get_FPY_OFFSET();

    /* Loop over this processor */
    for (iz = 0; iz < get_FPZ0_GRID(); iz++)
    {
        t1 = 0.0;
        for (iy = 0; iy < get_FPY0_GRID(); iy++)
           for (ix = 0; ix < get_FPX0_GRID(); ix++)
                if(ix == 0 && iy ==0 && pxoff == 0 && pyoff == 0)
                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[iz + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNZ_GRID();
    rmg::reduce(zvec.data(), ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        rmg::printlog("\n\n Planar average of the electrostatic density\n");
        rmg::printlog("&& %s\n", ab);
        for (ix = 0; ix < get_FNZ_GRID(); ix++)
        {
            t1 = ix * get_hzzgrid();
            //rmg::printlog(" %d %f %s\n", ix, zvec[ix] / get_FNX_GRID() / get_FNZ_GRID(), ab);
            rmg::printlog(" %d %e %s\n", ix, zvec[ix], ab);
        }
        rmg::printlog(" & %s\n", ab);
        fflush(NULL);
    }

}                               /* end get_avgd */

