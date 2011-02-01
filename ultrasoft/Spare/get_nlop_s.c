/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_nlop_s.c *****
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
 *   void get_nlop_s(ION *iptr, REAL *rtptr2, int ip, int icount, int *dvec)
 *   Generates the s-projectors 
 * INPUTS
 *   iptr: point to ion structure (see main.h)
 *   ip:   momentum l
 *   icount: number of grid point in the sphere
 *   dvec:  true indicates the point is in the sphere
 *          fault indicate not
 * OUTPUT
 *   rtptr2:  projectors
 * PARENTS
 *   get_nlop_smp.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

/* Generates s-projectors */
void get_nlop_s (ION * iptr, REAL * rtptr, int ip, int icount, int *dvec)
{

    int idx, ix, iy, iz, docount;
    REAL r, ax[3], xc, yc, zc;
    REAL invdr, t1;
    SPECIES *sp;


    /* Get species type */
    sp = &ct.sp[iptr->species];
    invdr = ONE / sp->drnlig;

    idx = docount = 0;
    xc = iptr->nlxcstart;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        yc = iptr->nlycstart;
        for (iy = 0; iy < sp->nldim; iy++)
        {

            zc = iptr->nlzcstart;
            for (iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {

                    ax[0] = xc - iptr->xtal[0];
                    ax[1] = yc - iptr->xtal[1];
                    ax[2] = zc - iptr->xtal[2];

                    r = metric (ax);
                    t1 = linint (&sp->betalig[ip][0], r, invdr);
                    rtptr[docount] = t1 * sqrt (1.0 / (4.0 * PI));

                    docount++;

                }               /* end if */

                idx++;

                zc += ct.hzgrid;
            }                   /* end for */

            yc += ct.hygrid;
        }                       /* end for */

        xc += ct.hxgrid;
    }                           /* end for */

    if (docount != icount)
        error_handler ("Problem with non-local generation");

}                               /* end get_nlop_s */

/******/
