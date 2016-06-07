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



/*
 * FUNCTION
 *   int claim_ion (int gridpe, ION * iptr,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid)
 *   This function determined which processor owns an ion. This depends on its spatial position relative
 *   to regions belonging to processors. Obviosuly, each ion has to belong to some processor and no two (or more)
 *   processors can share an ion. The function returns rank of the owner.
 *
 */


#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <math.h>


int claim_ion (double *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid)
{

    int  pe;
    int igridx, igridy, igridz;
    int xnode, ynode, znode;
    double t1, t2;



    /*Figure out grid coordinates of a grid point closest to position of ion
     * under consideration*/
    t1 = (xtal[0]) * (double) nxgrid;
    t1 = modf (t1, &t2);
    igridx = (int) t2;
    if (t1 > 0.5)
        igridx++;
    
    if (igridx >= nxgrid) igridx -= nxgrid ;
    
    t1 = (xtal[1]) * (double) nygrid;
    t1 = modf (t1, &t2);
    igridy = (int) t2;
    if (t1 > 0.5)
        igridy++;
    
    if (igridy >= nygrid) igridy -= nygrid ;
    
    t1 = (xtal[2]) * (double) nzgrid;
    t1 = modf (t1, &t2);
    igridz = (int) t2;
    if (t1 > 0.5)
        igridz++;
    
    if (igridz >= nzgrid) igridz -= nzgrid ;

    /*Now find the rank of the owner*/
//    pe = (igridx / pxgrid) * PE_Y * PE_Z + (igridy / pygrid) * PE_Z + (igridz / pzgrid);
    find_grid_owner(igridx, igridy, igridz, nxgrid, nygrid, nzgrid, &xnode, &ynode, &znode);

    pe = xnode * get_PE_Y() * get_PE_Z() +
         ynode * get_PE_Z() +
         znode;

    return (pe);

}                               /* end get_index */

/******/
