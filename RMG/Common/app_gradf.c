/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_gradf.c *****
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
 *   void app_gradf(FS0_GRID *f, FP0_GRID *wx, FP0_GRID *wy, FP0_GRID *wz)
 *   Generates the gradient of function defined on the FINE grid.
 * INPUTS
 *   S0_GRID *f: the function to be applied (see main.h)
 * OUTPUT
 *   P0_GRID *wx: (\partial/\partial x)f
 *   P0_GRID *wy: (\partial/\partial y)f
 *   P0_GRID *wz: (\partial/\partial z)f
 * PARENTS
 *   mg_eig_state.c subdiag_mpi.c subdiag_smp.c xcgga.c
 * CHILDREN
 *   trade_images2.c
 * SOURCE
 */


#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

void app_gradf (REAL * rho, FP0_GRID * wx, FP0_GRID * wy, FP0_GRID * wz)
{

    int iz, ix, iy;
    REAL t1, t2, t1x, t2x, t1y, t2y, t1z, t2z;
    REAL time1, time2, *rptr;
    REAL *wxr, *wyr, *wzr;
    FSS0_GRID *SSSptr;
    int ixs, iys;
    int ix1, iy1;

    ixs = (FPY0_GRID + 4) * (FPZ0_GRID + 4);
    iys = (FPZ0_GRID + 4);
    ix1 = FPY0_GRID * FPZ0_GRID;
    iy1 = FPZ0_GRID;
    wxr = (REAL *) wx;
    wyr = (REAL *) wy;
    wzr = (REAL *) wz;

    my_malloc (rptr, (FPX0_GRID + 4) * (FPY0_GRID + 4) * (FPZ0_GRID + 4), REAL);


    time1 = my_crtc ();
    /*trade_images2f(f, SSSptr); */
    /*trade_imagesx (f->s2, rptr, FPX0_GRID, FPY0_GRID, FPZ0_GRID, 2);*/
    trade_imagesx (rho, rptr, FPX0_GRID, FPY0_GRID, FPZ0_GRID, 2);

    switch (ct.ibrav)
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:
        t1x = 8.0 / (12.0 * ct.hxxgrid * ct.xside);
        t2x = 1.0 / (12.0 * ct.hxxgrid * ct.xside);

        t1y = 8.0 / (12.0 * ct.hyygrid * ct.yside);
        t2y = 1.0 / (12.0 * ct.hyygrid * ct.yside);

        t1z = 8.0 / (12.0 * ct.hzzgrid * ct.zside);
        t2z = 1.0 / (12.0 * ct.hzzgrid * ct.zside);

        for (ix = 2; ix < FPX0_GRID + 2; ix++)
        {

            for (iy = 2; iy < FPY0_GRID + 2; iy++)
            {

                for (iz = 2; iz < FPZ0_GRID + 2; iz++)
                {

                    wxr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2x * rptr[(ix - 2) * ixs + iy * iys + iz] +
                        -t1x * rptr[(ix - 1) * ixs + iy * iys + iz] +
                        t1x * rptr[(ix + 1) * ixs + iy * iys + iz] +
                        -t2x * rptr[(ix + 2) * ixs + iy * iys + iz];

                    wyr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2y * rptr[ix * ixs + (iy - 2) * iys + iz] +
                        -t1y * rptr[ix * ixs + (iy - 1) * iys + iz] +
                        t1y * rptr[ix * ixs + (iy + 1) * iys + iz] +
                        -t2y * rptr[ix * ixs + (iy + 2) * iys + iz];

                    wzr[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                        t2z * rptr[ix * ixs + iy * iys + iz - 2] +
                        -t1z * rptr[ix * ixs + iy * iys + iz - 1] +
                        t1z * rptr[ix * ixs + iy * iys + iz + 1] +
                        -t2z * rptr[ix * ixs + iy * iys + iz + 2];



                }               /* end for */
            }                   /* end for */
        }                       /* end for */

        break;

    case HEXAGONAL:

    
	/* This case may not be working, for some reason HEXAGONAL uses 3D array, while all 
	 * other cases use only 1D, this may be a problem if unions are not properly defined
	 * In either case, it would be best to have this case rewritten as 1D array to match 
	 * the above*/
	SSSptr = (FSS0_GRID *) rptr;
        
	t1 = sqrt (3.0) / 2.0;
        t1x = 8.0 / (12.0 * ct.hxxgrid * ct.xside);
        t2x = 1.0 / (12.0 * ct.hxxgrid * ct.xside);

        t1y = 8.0 / (t1 * 12.0 * ct.hxxgrid * ct.xside);
        t2y = 1.0 / (t1 * 12.0 * ct.hxxgrid * ct.xside);

        t1z = 8.0 / (12.0 * ct.hzzgrid * ct.zside);
        t2z = 1.0 / (12.0 * ct.hzzgrid * ct.zside);

        for (ix = 2; ix < FPX0_GRID + 2; ix++)
        {

            for (iy = 2; iy < FPY0_GRID + 2; iy++)
            {

                for (iz = 2; iz < FPZ0_GRID + 2; iz++)
                {

                    wx->s1.b[ix - 2][iy - 2][iz - 2] =
                        t2x * SSSptr->b[ix - 2][iy][iz] +
                        -t1x * SSSptr->b[ix - 1][iy][iz] +
                        t1x * SSSptr->b[ix + 1][iy][iz] + -t2x * SSSptr->b[ix + 2][iy][iz];

                    t1 = (SSSptr->b[ix + 1][iy - 1][iz] + SSSptr->b[ix][iy - 1][iz]) / 2.0;

                    t2 = (SSSptr->b[ix][iy + 1][iz] + SSSptr->b[ix - 1][iy + 1][iz]) / 2.0;

                    wy->s1.b[ix - 2][iy - 2][iz - 2] =
                        t2y * SSSptr->b[ix + 1][iy - 2][iz] +
                        -t1y * t1 + t1y * t2 + -t2y * SSSptr->b[ix - 1][iy + 2][iz];

                    wz->s1.b[ix - 2][iy - 2][iz - 2] =
                        t2z * SSSptr->b[ix][iy][iz - 2] +
                        -t1z * SSSptr->b[ix][iy][iz - 1] +
                        t1z * SSSptr->b[ix][iy][iz + 1] + -t2z * SSSptr->b[ix][iy][iz + 2];



                }               /* end for */
            }                   /* end for */
        }                       /* end for */
        break;

    default:
        error_handler ("Lattice type not implemented");
    }                           /* end switch */

    my_free (rptr);

    time2 = my_crtc ();
    rmg_timings (APPGRAD_TIME, (time2 - time1));



}                               /* end app_grad */


/******/
