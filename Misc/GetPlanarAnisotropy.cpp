/*
 *
 * Copyright 2022 The RMG Project Developers. See the COPYRIGHT file 
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
#include <iostream>
#include <fstream>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "transition.h"

double GetPlanarAnisotropy(void)
{
    double t[3];
    Rmg_L.cross_product(Rmg_L.a0, Rmg_L.a1, t);
    double xy_density = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NX_GRID(1)*Rmg_G->get_NY_GRID(1));
    xy_density = 1.0 / xy_density;

    Rmg_L.cross_product(Rmg_L.a0, Rmg_L.a2, t);
    double xz_density = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NX_GRID(1)*Rmg_G->get_NZ_GRID(1));
    xz_density = 1.0 / xz_density;

    Rmg_L.cross_product(Rmg_L.a1, Rmg_L.a2, t);
    double yz_density = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NY_GRID(1)*Rmg_G->get_NZ_GRID(1));
    yz_density = 1.0 / yz_density;

    if(ct.verbose && pct.gridpe == 0) 
        rmg_printf("Planar Density = %12.6f  %12.6f  %12.6f\n", xy_density, xz_density, yz_density);
    
    double dmax = std::max(xy_density, xz_density);
    dmax = std::max(dmax, yz_density);
    double dmin = std::min(xy_density, xz_density);
    dmin = std::min(dmin, yz_density);
    if(ct.verbose && pct.gridpe == 0) 
        rmg_printf("Planar Anisotropy = %12.6f\n", dmax/dmin);

    return dmax / dmin;
}
