/****f* QMD-MGDFT/init_efield.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2001  Vincent Meunier, Wenchang Lu, Jerzy Bernholc
 * FUNCTION
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

#include "main.h"

#define SMALL_EFIELD 1.0e-20

void init_efield (double * vnuc, double efield[3])
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
	    rmg_printf ("\n EXTERNAL ELECTRIC FIELD INITIALIZATION ");
	    rmg_printf ("\n ====================================== \n");
	    rmg_printf ("\n Electric Field: (%12.8f, %12.8f, %12.8f) Ha/bohr\n",
		    efield[0], efield[1], efield[2]);
    /* find the cartesian coordinate of the corner */
    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * get_hxgrid() * get_xside() * get_PX0_GRID();
    yoff = iy * get_hygrid() * get_yside() * get_PY0_GRID();
    zoff = iz * get_hzgrid() * get_zside() * get_PZ0_GRID();
    incix = get_FPY0_GRID() * get_FPZ0_GRID();
    inciy = get_FPZ0_GRID();
/*   xoff, yoff, zoff are same for fine and coarse grid  */

/*  we need to add the saw-tooth potential in fine grid */
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
                vnuc[idx] += rx * efield[0] + ry * efield[1] + rz * efield[2];
            }
        }
    }

}                               /* end init_efield */

/******/
