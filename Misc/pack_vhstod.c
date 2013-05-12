/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/pack_vhstod.c *****
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
 *   void pack_vhstod(rmg_double_t *s, REAL *d, int dimx, int dimy, int dimz)
 *   Used to pack grids when computing the hartree potential 
 *   For periodic system, d = s
 *   For surface system,  d[x,y,Nz/2+1:Nz/2+Nz] = s[x,y,1:Nz]
 *   FOr cluster system:  d[Nx/2+1:Nx/2+Nx, ...] = s[1:Nx, 1:Ny, 1:Nz] 
 * INPUTS
 *   s[dmix*dimy*dimz]
 *   dimx, dimy, dimz: dimension of the wave function grid
 * OUTPUT
 *   d[ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid]
 *       see init.c for ct.vh_pxgrid, ...
 *       array in the hatree potential grid
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include "main.h"
#include <float.h>
#include <math.h>

/* This function is used to pack grids when computing the hartree potential */
void pack_vhstod (rmg_double_t * s, REAL * d, int dimx, int dimy, int dimz)
{
    int ix, iy, iz;
    int pex, pey, pez;
    int sxlo, sylo, szlo;
    int sxhi, syhi, szhi;
    int dxlo, dylo, dzlo;
    int dxhi, dyhi, dzhi;
    int loidx, hiidx;
    int hix, hiy, hiz;
    int hincx, hincy;

    Dprintf ("Checking boundaryflag value of %d", ct.boundaryflag);
    if (ct.boundaryflag == PERIODIC)
    {

        for (ix = 0; ix < (dimx * dimy * dimz); ix++)
            d[ix] = s[ix];
        return;

    }                           /* end if */

    pe2xyz (pct.gridpe, &pex, &pey, &pez);
    sxlo = pex * dimx;
    sylo = pey * dimy;
    szlo = pez * dimz;
    sxhi = sxlo + dimx - 1;
    syhi = sylo + dimy - 1;
    szhi = szlo + dimz - 1;

    dxlo = pex * 2 * dimx - ct.psi_nxgrid / 2;
    dylo = pey * 2 * dimy - ct.psi_nygrid / 2;
    dzlo = pez * 2 * dimz - ct.psi_nzgrid / 2;
    dxhi = dxlo + 2 * dimx - 1;
    dyhi = dylo + 2 * dimy - 1;
    dzhi = dzlo + 2 * dimz - 1;
    hincy = 2 * dimz;
    hincx = 4 * dimy * dimz;

    if (ct.boundaryflag == SURFACE)
    {

        dxlo = sxlo;
        dylo = sylo;
        dxhi = sxhi;
        dyhi = syhi;
        hincx = 2 * dimy * dimz;

    }                           /* end if */


    loidx = 0;
    for (ix = sxlo; ix <= sxhi; ix++)
    {

        for (iy = sylo; iy <= syhi; iy++)
        {

            for (iz = szlo; iz <= szhi; iz++)
            {


                /* Check if we map into the double grid */
                if ((ix >= dxlo) && (ix <= dxhi) &&
                    (iy >= dylo) && (iy <= dyhi) && (iz >= dzlo) && (iz <= dzhi))
                {

                    hix = ix - dxlo;
                    hiy = iy - dylo;
                    hiz = iz - dzlo;
                    hiidx = hix * hincx + hiy * hincy + hiz;

                    /* Yes there is a mapping so transfer the data */
                    d[hiidx] = s[loidx];

                }               /* end if */

                loidx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_vhstod */

/******/
