/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_index.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   int claim_ion (int gridpe, ION * iptr,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid)
 *   This function determined which processor owns an ion. This depends on its spatial position relative
 *   to regions belonging to processors. Obviosuly, each ion has to belong to some processor and no two (or more)
 *   processors can share an ion. The function returns rank of the owner.
 *
 * SOURCE
 */


#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <math.h>


int claim_ion (rmg_double_t *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid)
{

    int  pe;
    int ii, jj, kk, ilow, ihi, klow, khi, jlow, jhi;
    int igridx, igridy, igridz;
    int xnode, ynode, znode;
    rmg_double_t t1, t2;



    /*Figure out grid coordinates of a grid point closest to position of ion
     * under consideration*/
    t1 = (xtal[0]) * (rmg_double_t) nxgrid;
    t1 = modf (t1, &t2);
    igridx = (int) t2;
    if (t1 > 0.5)
        igridx++;
    
    if (igridx >= nxgrid) igridx -= nxgrid ;
    
    t1 = (xtal[1]) * (rmg_double_t) nygrid;
    t1 = modf (t1, &t2);
    igridy = (int) t2;
    if (t1 > 0.5)
        igridy++;
    
    if (igridy >= nygrid) igridy -= nygrid ;
    
    t1 = (xtal[2]) * (rmg_double_t) nzgrid;
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
