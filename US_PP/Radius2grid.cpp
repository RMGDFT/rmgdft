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
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Atomic.h"


/*For a quantity localized around ionic positions, this function finds radius in number of grid points
 * given a radius in a.u.*/
int Radius2grid (double radius, double mingrid_spacing)
{

    double scale, t1, t2;
    int it1, dim, ibrav;

    ibrav = Rmg_L.get_ibrav_type();
    
    /* Set the scaling factor for determining the radius of the local grids */
    scale = 1.0;
    if (ibrav == CUBIC_BC)
	scale = 1.1;
    if (ibrav == CUBIC_FC)
	scale = 1.3;
        
    t1 = 2.0 * scale * radius / mingrid_spacing;
    t1 = modf (t1, &t2);
    it1 = (int) t2;
    if (t1 > 0.5)
          it1++;
    if (!(it1 % 2))
        it1++;
    dim = it1;

    return dim;
}
