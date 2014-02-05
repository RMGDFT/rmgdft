/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_phase.c *****
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
 *   void get_phase(ION *iptr, rmg_double_t *rtptr, int ip, int icount, int *dvec)
 *   Generates the phase factors for the k-points 
 * INPUTS
 *   iptr: point to structure ION (see main.h )
 *   ip:   momentum l
 *   icount: number of grid point in the sphere
 *   dvec:  true indicates the point is in the sphere
 *          fault indicate not
 * OUTPUT
 *   rtptr: pct.phaseptr (see main.h )
 * PARENTS
 *   get_nloc_smp.c
 * CHILDREN
 *   to_cartesian.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "common_prototypes.h"

void get_phase (ION * iptr, rmg_double_t * rtptr, int ip, int icount, int *dvec)
{

    int kpt, idx, ix, iy, iz, docount;
    rmg_double_t ax[3], bx[3], xc, yc, zc;
    rmg_double_t kdr;
    rmg_double_t hxgrid, hygrid, hzgrid;
    SPECIES *sp;

    hxgrid = get_hxgrid();
    hygrid = get_hygrid();
    hzgrid = get_hzgrid();

    if (rtptr == NULL)
        return;

    /* Get species type */
    sp = &ct.sp[iptr->species];


    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

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
                        kdr = ct.kp[kpt].kvec[0] * bx[0] +
                            ct.kp[kpt].kvec[1] * bx[1] + ct.kp[kpt].kvec[2] * bx[2];


                        rtptr[docount + 2 * kpt * icount] = cos (kdr);
                        rtptr[docount + 2 * kpt * icount + icount] = sin (kdr);
                        docount++;

                    }           /* end if */

                    idx++;


                    zc += hzgrid;
                }               /* end for */

                yc += hygrid;
            }                   /* end for */

            xc += hxgrid;
        }                       /* end for */

        if (docount != icount)
            error_handler ("Something wrong with docount");
    }                           /* end for */


}                               /* end get_phase */

/******/
