/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_nlop_d.c *****
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
 *   void get_nlop_d(ION *iptr, REAL *rtptr2, int ip, int icount, int *dvec)
 *   Generates the d-projectors 
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

void get_nlop_d (ION * iptr, REAL * rtptr, int ip, int icount, int *dvec)
{

    int idx, ix, iy, iz, docount;
    REAL r, ax[3], bx[3], xc, yc, zc, t1, cc, t2, rsq1;
    REAL *r1, *r2, *r3, *r4, *r5, invdr;
    SPECIES *sp;

    r1 = rtptr;
    r2 = r1 + icount;
    r3 = r2 + icount;
    r4 = r3 + icount;
    r5 = r4 + icount;

    /* Get species type */
    sp = &ct.sp[iptr->species];
    invdr = ONE / sp->drnlig;

    cc = sqrt (5.0 / (4.0 * PI));
    t2 = sqrt (3.0);

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
                    to_cartesian (ax, bx);
                    r = metric (ax);
                    rsq1 = r * r + 1.0e-20;
                    t1 = linint (&sp->betalig[ip][0], r, invdr);
                    t1 = t1 * t2;
                    r1[docount] = cc * t1 * bx[0] * bx[1] / rsq1;
                    r2[docount] = cc * t1 * bx[0] * bx[2] / rsq1;
                    r3[docount] = cc * t1 * (t2 * bx[2] * bx[2] - rsq1 / t2) / (2.0 * rsq1);
                    r4[docount] = cc * t1 * bx[1] * bx[2] / rsq1;
                    r5[docount] = cc * t1 * (bx[0] * bx[0] - bx[1] * bx[1]) / (2.0 * rsq1);

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

}                               /* end get_nlop_d */

/******/
