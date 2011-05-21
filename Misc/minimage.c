/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/minimage.c *****
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
 *   REAL minimage(ION *ip1, ION *ip2, REAL *xtal_r)
 *   This function returns the minimum distance between two ions in the
 *   supercell. It also returns the x, y, and z coordinates for the minimum
 *   distance. 
 * INPUTS
 *   ip1, ip2: point to structure ION, only use the coordinates
 * OUTPUT
 *   return the minimum distance
 *   xtal_r: store the shortest vector 
 * PARENTS
 *   get _te.c iiforce.c
 * CHILDREN
 *   metric.c
 * SOURCE
 */





#include "main.h"
#include <float.h>
#include <math.h>


REAL minimage (ION * ip1, ION * ip2, REAL * xtal_r)
{

    int ix, iy, iz, idx, idxmin = 0;
    REAL r[27], ax[3], rmin, x[27], y[27], z[27];



    /* If cluster boundary conditions just compute r and return */
    if (ct.boundaryflag == CLUSTER)
    {

        xtal_r[0] = ip1->xtal[0] - ip2->xtal[0];
        xtal_r[1] = ip1->xtal[1] - ip2->xtal[1];
        xtal_r[2] = ip1->xtal[2] - ip2->xtal[2];

        return metric (xtal_r);

    }                           /* end if */


    /* Loop over all possible combinations and generate r */

    idx = 0;
    for (ix = -1; ix <= 1; ix++)
    {

        for (iy = -1; iy <= 1; iy++)
        {

            for (iz = -1; iz <= 1; iz++)
            {

                x[idx] = ip1->xtal[0] - ip2->xtal[0] + (REAL) ix;
                y[idx] = ip1->xtal[1] - ip2->xtal[1] + (REAL) iy;
                z[idx] = ip1->xtal[2] - ip2->xtal[2] + (REAL) iz;
                if (ct.boundaryflag == SURFACE)
                {

                    z[idx] = ip1->xtal[2] - ip2->xtal[2];

                }               /* end if */

                ax[0] = x[idx];
                ax[1] = y[idx];
                ax[2] = z[idx];

                r[idx] = metric (ax);

                idx++;

            }                   /* end for ix */

        }                       /* end for iy */

    }                           /* end for iz */


    /* Next we find rmin */
    rmin = 100000.0;
    for (idx = 0; idx < 27; idx++)
    {

        /* Gaurd against roundoff */
        if ((r[idx] > 1.0e-5) && (r[idx] < rmin))
        {

            rmin = r[idx];
            idxmin = idx;

        }                       /* end if */

    }                           /* end for */


    xtal_r[0] = x[idxmin];
    xtal_r[1] = y[idxmin];
    xtal_r[2] = z[idxmin];

    return rmin;

}                               /* end minimage */

/******/
