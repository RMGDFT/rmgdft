/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_smooth.c *****
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
 *   void app_smooth(S0_GRID *f, S0_GRID *work, REAL sfac)
 *   Smooths the function f by averaging the neighbors 
 * INPUTS
 *   f: array to be smoothed
 *   sfac: scale factor, it is useless now
 * OUTPUT
 *   work: smoothed results 
 * PARENTS
 *   mg_eig_state.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>


/**

 Smooths f and returns result in work
*/
void app_smooth (S0_GRID * f, S0_GRID * work, REAL sfac)
{

    int iz, ix, iy;
    REAL scale, ec, fc, crn, cc, temp;

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;


    for (ix = 1; ix <= PX0_GRID; ix++)
    {

        for (iy = 1; iy <= PY0_GRID; iy++)
        {

            for (iz = 1; iz <= PZ0_GRID; iz++)
            {


                temp = cc * f->s1.b[ix][iy][iz] +
                    fc * (f->s1.b[ix - 1][iy][iz] +
                          f->s1.b[ix + 1][iy][iz] +
                          f->s1.b[ix][iy - 1][iz] +
                          f->s1.b[ix][iy + 1][iz] +
                          f->s1.b[ix][iy][iz - 1] + f->s1.b[ix][iy][iz + 1]);

                temp +=
                    ec * (f->s1.b[ix - 1][iy][iz - 1] +
                          f->s1.b[ix + 1][iy][iz - 1] +
                          f->s1.b[ix][iy - 1][iz - 1] +
                          f->s1.b[ix][iy + 1][iz - 1] +
                          f->s1.b[ix - 1][iy - 1][iz] +
                          f->s1.b[ix - 1][iy + 1][iz] +
                          f->s1.b[ix + 1][iy - 1][iz] +
                          f->s1.b[ix + 1][iy + 1][iz] +
                          f->s1.b[ix - 1][iy][iz + 1] +
                          f->s1.b[ix + 1][iy][iz + 1] +
                          f->s1.b[ix][iy - 1][iz + 1] + f->s1.b[ix][iy + 1][iz + 1]);

                temp +=
                    crn * (f->s1.b[ix - 1][iy - 1][iz - 1] +
                           f->s1.b[ix - 1][iy - 1][iz + 1] +
                           f->s1.b[ix - 1][iy + 1][iz - 1] +
                           f->s1.b[ix - 1][iy + 1][iz + 1] +
                           f->s1.b[ix + 1][iy - 1][iz - 1] +
                           f->s1.b[ix + 1][iy - 1][iz + 1] +
                           f->s1.b[ix + 1][iy + 1][iz - 1] + f->s1.b[ix + 1][iy + 1][iz + 1]);

                work->s1.b[ix][iy][iz] = scale * temp;


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



}                               /* end app_smooth */






/******/
