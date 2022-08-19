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

double GetPlanarAnisotropy(double *density)
{
    double t[3];
    Rmg_L.cross_product(Rmg_L.a0, Rmg_L.a1, t);
    density[0] = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NX_GRID(1)*Rmg_G->get_NY_GRID(1));
    density[0] = 1.0 / density[0];

    Rmg_L.cross_product(Rmg_L.a0, Rmg_L.a2, t);
    density[1] = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NX_GRID(1)*Rmg_G->get_NZ_GRID(1));
    density[1] = 1.0 / density[1];

    Rmg_L.cross_product(Rmg_L.a1, Rmg_L.a2, t);
    density[2] = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) / 
	                (Rmg_G->get_NY_GRID(1)*Rmg_G->get_NZ_GRID(1));
    density[2] = 1.0 / density[2];

    double dmax = std::max(density[0], density[1]);
    dmax = std::max(dmax, density[2]);
    double dmin = std::min(density[0], density[1]);
    dmin = std::min(dmin, density[2]);

    return dmax / dmin;
}
