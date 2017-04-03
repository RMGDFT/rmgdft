
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
#include <string.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgException.h"
#include "Functional.h"
#include "Solvers.h"
#include "RmgParallelFft.h"

void GetPhaseSpecies (SPECIES *sp , std::complex<double> *rtptr)
{

    int kpt, idx, ix, iy, iz;
    int dimx, dimy, dimz, pbasis;
    double ax[3], bx[3], xc, yc, zc;
    double kdr;
    double hxgrid, hygrid, hzgrid;


    if (rtptr == NULL)
        return;

    hxgrid = Rmg_G->get_hxgrid(1);
    hygrid = Rmg_G->get_hygrid(1);
    hzgrid = Rmg_G->get_hzgrid(1);


    dimx = sp->nldim;
    dimy = sp->nldim;
    dimz = sp->nldim;

    pbasis = dimx * dimy * dimz;


    for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        int kpt1 = kpt + pct.kstart;
        for (ix = 0; ix < dimx; ix++)
        {

            xc = (ix-dimx/2) * hxgrid;

            for (iy = 0; iy < dimy; iy++)
            {
                yc = (iy-dimy/2) * hygrid;

                for (iz = 0; iz < dimz; iz++)
                {
                    zc = (iz-dimz/2) * hzgrid;

                    ax[0] = xc ;
                    ax[1] = yc ;
                    ax[2] = zc ;
                    to_cartesian (ax, bx);
                    kdr = ct.kp[kpt1].kvec[0] * bx[0] +
                        ct.kp[kpt1].kvec[1] * bx[1] + ct.kp[kpt1].kvec[2] * bx[2];

                    idx = ix * dimy * dimz + iy * dimz + iz;

                    rtptr[idx + kpt * pbasis] = exp(std::complex<double>(0.0, -kdr));


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end get_phase */

