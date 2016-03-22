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
 *   void init_efield(double *vnuc)
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
#include "common_prototypes.h"
#include "grid.h"
#include "main.h"

#define SMALL_EFIELD 1.0e-20

void init_efield (double * vnuc)
{

    double norm_field_0;
    int ix, iy, iz;
    int jx, jy, jz;
    int idx;
    int incix, inciy;
    double xoff, yoff, zoff;
    double rx, ry, rz;
    double d_field;
    double hxxgrid, hyygrid, hzzgrid;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();


    /* normalization of the electric field direction */
    norm_field_0 = sqrt (pow (ct.x_field_0 + SMALL_EFIELD, 2) + pow (ct.y_field_0 + SMALL_EFIELD, 2) + pow (ct.z_field_0 + SMALL_EFIELD, 2));
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
    xoff = ix * get_hxgrid() * get_xside() * get_PX0_GRID();
    yoff = iy * get_hygrid() * get_yside() * get_PY0_GRID();
    zoff = iz * get_hzgrid() * get_zside() * get_PZ0_GRID();
    incix = get_FPY0_GRID() * get_FPZ0_GRID();
    inciy = get_FPZ0_GRID();
/*   xoff, yoff, zoff are same for fine and coarse grid  */

/*  we need to add the saw-tooth potential in fine grid */
    if(ct.runflag == 5 || ct.runflag == 6) return;
    for (jx = 0; jx < get_FPX0_GRID(); jx++)
    {
        rx = jx * hxxgrid * get_xside() + xoff;
        for (jy = 0; jy < get_FPY0_GRID(); jy++)
        {
            ry = jy * hyygrid * get_yside() + yoff;
            for (jz = 0; jz < get_FPZ0_GRID(); jz++)
            {
                idx = jx * incix + jy * inciy + jz;
                rz = jz * hzzgrid * get_zside() + zoff;
                d_field = rx * ct.x_field_0 + ry * ct.y_field_0 + rz * ct.z_field_0;
                vnuc[idx] += ct.e_field * d_field;
            }
        }
    }

}                               /* end init_efield */

/******/
