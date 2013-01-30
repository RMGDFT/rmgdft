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
 *   int get_index(ION *iptr, int *Aix, int *Aiy, int *Aiz,
 *                 int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
 *                 int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
 *                 int nxgrid, int nygrid, int nzgrid,
 *                 REAL *xcstart, REAL *ycstart, REAL *zcstart)
 *   This function generates the indices needed to map a spherical 
 *   short-ranged operator onto a rectangular grid. This is required
 *   for the non-local potential operators. Since the code uses
 *   domain-decomposition with each individual processor mapped to
 *   a specific region of space the function has to take that into 
 *   account as well as taking care of periodic boundary conditions.
 * INPUTS
 *   iptr: Pointer to the ion control structure
 *   pxgrid: Number of points in the processor grid in the x-direction
 *   pygrid: Number of points in the processor grid in the y-direction
 *   pzgrid: Number of points in the processor grid in the z-direction
 *   nxgrid: Number of points in the global grid in the x-direction
 *   nygrid: Number of points in the global grid in the y-direction
 *   nzgrid: Number of points in the global grid in the z-direction
 *   xcstart: Physical x-coordinate of the lower corner of the processor grid
 *   ycstart: Physical y-coordinate of the lower corner of the processor grid
 *   zcstart: Physical z-coordinate of the lower corner of the processor grid
 * OUTPUT 
 *    returned value
 *   If true indicates that the spherical short-ranged operator
 *   actually exists in this particular processors region of space. If
 *   false indicates that the operator does not exist in this region
 *   of space.
 * PARENTS
 *   get_nlop_smp.c init_nuc.c lforce.c nlccforce.c nlforce.c
 * CHILDREN
 *   pe2xyz.c
 * SOURCE
 */


#include "main.h"
#include <math.h>


int test_overlap (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid)
{

    int idx, ix, iy, iz, ic, map, ilow, ihi, jlow, jhi, klow, khi;
    int ii, jj, kk, icenter, itmp;
    int ixstart, iystart, izstart, icut;
    int pxoffset, pyoffset, pzoffset;
    REAL t1, t2, vect [3];
    REAL xcstart, ycstart, zcstart;



    /* Generate range of indices over which the operator */
    /* will be mapped onto the global grid.              */
    t1 = iptr->xtal[0] * (REAL) nxgrid;
    t1 = modf (t1, &t2);
    ic = (int) t2;
    if (t1 > 0.5)
        ic++;
    ixstart = ic - cdim / 2;
    xcstart = ixstart / (REAL) nxgrid;
    


    t1 = iptr->xtal[1]  * (REAL) nygrid;
    t1 = modf (t1, &t2);
    ic = (int) t2;
    if (t1 > 0.5)
        ic++;
    iystart = ic - cdim / 2;
    ycstart = iystart / (REAL) nygrid;


    t1 = iptr->xtal[2] * (REAL) nzgrid;
    t1 = modf (t1, &t2);
    ic = (int) t2;
    if (t1 > 0.5)
        ic++;
    izstart = ic - cdim / 2;
    zcstart = izstart / (REAL) nzgrid;



    /* Generate range of indices over which the projectors */
    /* will be mapped onto the global grid.                */
    idx = ixstart;
    for (ix = 0; ix < cdim; ix++)
    {

        if (idx < 0)
        {

            Aix[ix] = nxgrid + idx;

        }
        else if (idx == 0)
        {

            Aix[ix] = idx;

        }
        else if (idx > (nxgrid - 1))
        {

            Aix[ix] = idx - nxgrid;

        }
        else
        {

            Aix[ix] = idx;

        }                       /* end if */

        idx++;

    }                           /* end for */


    idx = iystart;
    for (iy = 0; iy < cdim; iy++)
    {

        if (idx < 0)
        {

            Aiy[iy] = nygrid + idx;

        }
        else if (idx == 0)
        {

            Aiy[iy] = 0;

        }
        else if (idx > (nygrid - 1))
        {

            Aiy[iy] = idx - nygrid;

        }
        else
        {

            Aiy[iy] = idx;

        }                       /* end if */

        idx++;

    }                           /* end for */


    idx = izstart;
    for (iz = 0; iz < cdim; iz++)
    {

        if (idx < 0)
        {

            Aiz[iz] = nzgrid + idx;

        }
        else if (idx == 0)
        {

            Aiz[iz] = idx;

        }
        else if (idx > (nzgrid - 1))
        {

            Aiz[iz] = idx - nzgrid;

        }
        else
        {

            Aiz[iz] = idx;

        }                       /* end if */

        idx++;

    }                           /* end for */


    /* Now we need to determine if any of this ions */
    /* projector maps onto this processors space.   */
    pe2xyz (gridpe, &ii, &jj, &kk);
    find_node_offsets(gridpe, nxgrid, nygrid, nzgrid,
                      &pxoffset, &pyoffset, &pzoffset);

    ilow = pxoffset;
    jlow = pyoffset;
    klow = pzoffset;
    ihi = ilow + pxgrid - 1;
    jhi = jlow + pygrid - 1;
    khi = klow + pzgrid - 1;

    ii = jj = kk = FALSE;
    for (idx = 0; idx < cdim; idx++)
    {

        if ((Aix[idx] >= ilow) && (Aix[idx] <= ihi))
            ii = TRUE;
        if ((Aiy[idx] >= jlow) && (Aiy[idx] <= jhi))
            jj = TRUE;
        if ((Aiz[idx] >= klow) && (Aiz[idx] <= khi))
            kk = TRUE;

    }                           /* end for */

    map = ii & jj;
    map = map & kk;

    if (!map)
        return map;
    
    /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
    /*xtal vector between ion and left bottom corner of the box */
    /*Then, substract vector between left bottom corner of the box and center of the box */
    vect[0] = iptr->xtal[0] -  xcstart;
    vect[0] -= (cdim / 2) / (REAL) nxgrid;
    
    vect[1] = iptr->xtal[1] -  ycstart;
    vect[1] -= (cdim / 2) / (REAL) nygrid;
    
    vect[2] = iptr->xtal[2] -  zcstart;
    vect[2] -= (cdim / 2) / (REAL) nzgrid;
        
        
    icenter = cdim / 2;
    icut = (icenter + 1) * (icenter + 1);
            
    map = 0;
    for (ix = 0; ix < cdim; ix++)
    {
        for (iy = 0; iy < cdim; iy++)
        {
            for (iz = 0; iz < cdim; iz++)
            {
                        
                if ((((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                            ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                            ((Aiz[iz] >= klow) && (Aiz[iz] <= khi))))
                {
                    /* Cut it off if required */
                    itmp =
                        (ix - icenter) * (ix - icenter) +
                        (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                    if (icut >= itmp) 
                    {
                        map++;
                        return (map);
                    }

                }
            }
        }
    }

    return (map);



}                               /* end get_index */

/******/
