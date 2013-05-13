/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/pack_vhdtos.c *****
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
 *   void pack_vhdtos(rmg_double_t *s, rmg_double_t *d, int dimx, int dimy, int dimz)
 *   Used to pack grids when computing the hartree potential 
 *   For periodic system, s = d
 *   For surface system, s[x,y,1:Nz] = d[x,y,Nz/2+1:Nz/2+Nz)
 *   FOr cluster system: s[1:Nx, 1:Ny, 1:Nz] = d[Nx/2+1:Nx/2+Nx, ...]
 * INPUTS
 *   d[ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid]
 *     see init.c for ct.vh_pxgrid, ...
 *   dimx, dimy, dimz: dimension of the wave function grid
 * OUTPUT
 *   s[dimx*dimy*dimz] array in the wave function grid
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */
#include "main.h"
#include <float.h>
#include <math.h>

void pack_vhdtos (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz)
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
            s[ix] = d[ix];
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
                    s[loidx] = d[hiidx];

                }               /* end if */

                loidx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_vhdtos */

/******/
