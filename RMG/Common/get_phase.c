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


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "common_prototypes.h"

void get_phase (ION * iptr, double * rtptr, int icount, int *dvec)
{

    int kpt, idx, ix, iy, iz, docount;
    int dimx, dimy, dimz, pbasis;
    double ax[3], bx[3], xc, yc, zc;
    double kdr;
    double hxgrid, hygrid, hzgrid;
    SPECIES *sp;

    hxgrid = get_hxgrid();
    hygrid = get_hygrid();
    hzgrid = get_hzgrid();

    if (rtptr == NULL)
        return;

    /* Get species type */
    sp = &ct.sp[iptr->species];

    dimx = sp->nldim;
    dimy = sp->nldim;
    dimz = sp->nldim;
    pbasis = dimx * dimy * dimz;


    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        idx = 0;
        for (ix = 0; ix < dimx; ix++)
        {
    
           xc = iptr->nlxcstart + ix * hxgrid;

            for (iy = 0; iy < dimy; iy++)
            {
                yc = iptr->nlycstart + iy * hygrid;

                for (iz = 0; iz < dimz; iz++)
                {
                   zc = iptr->nlzcstart + iz * hzgrid;

                   // ax[0] = xc - iptr->xtal[0];
                   // ax[1] = yc - iptr->xtal[1];
                   // ax[2] = zc - iptr->xtal[2];
                    ax[0] = xc ;
                    ax[1] = yc ;
                    ax[2] = zc ;
                    to_cartesian (ax, bx);
                    kdr = ct.kp[kpt].kvec[0] * bx[0] +
                          ct.kp[kpt].kvec[1] * bx[1] + ct.kp[kpt].kvec[2] * bx[2];

                    idx = ix * dimy * dimz + iy * dimz + iz;

                    rtptr[idx + 2 * kpt * pbasis] = cos (kdr);
                    rtptr[idx + 2 * kpt * pbasis + pbasis] = sin (kdr);


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end get_phase */

/******/
