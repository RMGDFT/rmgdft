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
 *   processors can share an ion
 *
 * SOURCE
 */


#include "main.h"
#include <math.h>


int claim_ion (int gridpe, ION * iptr,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid)
{

    int  map;
    int ii, jj, kk, ilow, ihi, klow, khi, jlow, jhi;
    int igridx, igridy, igridz;
    REAL t1, t2;



    /*Figure out grid coordinates of a grid point closest to position of ion
     * under consideration*/
    t1 = (iptr->xtal[0]) * (REAL) nxgrid;
    t1 = modf (t1, &t2);
    igridx = (int) t2;
    if (t1 > 0.5)
        igridx++;
    
    if (igridx >= nxgrid) igridx -= nxgrid ;
    
    t1 = (iptr->xtal[1]) * (REAL) nygrid;
    t1 = modf (t1, &t2);
    igridy = (int) t2;
    if (t1 > 0.5)
        igridy++;
    
    if (igridy >= nygrid) igridy -= nygrid ;
    
    t1 = (iptr->xtal[2]) * (REAL) nzgrid;
    t1 = modf (t1, &t2);
    igridz = (int) t2;
    if (t1 > 0.5)
        igridz++;
    
    if (igridz >= nzgrid) igridz -= nzgrid ;



    /* Now we need to determine whether the grid point belong to
    /* belongs to current processor  */
    pe2xyz (gridpe, &ii, &jj, &kk);
    ilow = ii * pxgrid;
    jlow = jj * pygrid;
    klow = kk * pzgrid;
    ihi = ilow + pxgrid - 1;
    jhi = jlow + pygrid - 1;
    khi = klow + pzgrid - 1;


    ii = jj = kk = FALSE;

    if ((igridx >= ilow) && (igridx <= ihi))
	ii = TRUE;
    if ((igridy >= jlow) && (igridy <= jhi))
	jj = TRUE;
    if ((igridz >= klow) && (igridz <= khi))
	kk = TRUE;


    map = ii & jj;
    map = map & kk;
    return map;

}                               /* end get_index */

/******/
