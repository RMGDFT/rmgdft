/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/init_efield.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2001  Vincent Meunier, Wenchang Lu, Jerzy Bernholc
 * FUNCTION
 *   void init_efield(REAL *vnuc)
 *   Set up the electric field (add it to vnuc already defined)
 *   sawtooth potential
 * INPUTS
 *   vnuc: local part of pseudopotential
 * OUTPUT
 *   vnuc: local part of pseudopotential + V_efield 
 * PARENTS
 *   init_nuc.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void init_efield (REAL * vnuc)
{

    REAL norm_field_0;
    int ix, iy, iz;
    int jx, jy, jz;
    int idx;
    int incix, inciy;
    REAL xoff, yoff, zoff;
    REAL rx, ry, rz;
    REAL d_field;

    /* normalization of the electric field direction */
    norm_field_0 = sqrt (pow (ct.x_field_0, 2) + pow (ct.y_field_0, 2) + pow (ct.z_field_0, 2));
    ct.x_field_0 = ct.x_field_0 / norm_field_0;
    ct.y_field_0 = ct.y_field_0 / norm_field_0;
    ct.z_field_0 = ct.z_field_0 / norm_field_0;
    if (pct.gridpe == 0)
    {
        printf ("\n EXTERNAL ELECTRIC FIELD INITIALIZATION ");
        printf ("\n ====================================== \n");
        printf ("\n E_field = %16.8f", ct.e_field);
        printf ("\n Electric Field Direction: (%12.8f, %12.8f, %12.8f)",
                ct.x_field_0, ct.y_field_0, ct.z_field_0);
        printf ("\n Potential Slope: %16.8f Hartree/a.u", ct.e_field);
        printf ("\n Potential Slope: %16.8f eV/A\n", ct.e_field * Ha_eV / a0_A);
    }
    /* find the cartesian coordinate of the corner */
    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * ct.hxgrid * ct.xside * pct.PX0_GRID;
    yoff = iy * ct.hygrid * ct.yside * pct.PY0_GRID;
    zoff = iz * ct.hzgrid * ct.zside * pct.PZ0_GRID;
    incix = pct.FPY0_GRID * pct.FPZ0_GRID;
    inciy = pct.FPZ0_GRID;
/*   xoff, yoff, zoff are same for fine and coarse grid  */

/*  we need to add the saw-tooth potential in fine grid */
    for (jx = 0; jx < pct.FPX0_GRID; jx++)
    {
        rx = jx * ct.hxxgrid * ct.xside + xoff;
        for (jy = 0; jy < pct.FPY0_GRID; jy++)
        {
            ry = jy * ct.hyygrid * ct.yside + yoff;
            for (jz = 0; jz < pct.FPZ0_GRID; jz++)
            {
                idx = jx * incix + jy * inciy + jz;
                rz = jz * ct.hzzgrid * ct.zside + zoff;
                d_field = rx * ct.x_field_0 + ry * ct.y_field_0 + rz * ct.z_field_0;
                vnuc[idx] += ct.e_field * d_field;
            }
        }
    }

}                               /* end init_efield */

/******/
