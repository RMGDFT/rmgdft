/************************** SVN Revision Information **************************
 **    $Id: linint.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/linint.c *****
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
 *   REAL linint(REAL *y, REAL r, REAL invdr)
 *   Interpolation routine for data defined on a linear grid.
 * INPUTS
 *   y:     array of y-values
 *   rv:    r value being interpolated to
 *   invdr: inverse linear grid spacing
 * OUTPUT
 *   interpolated value is returned
 * PARENTS
 *   get_nlop_d.c get_nlop_p.c get_nlop_s.c init_nuc.c init_wflcao.c
 *   lforce.c nlccforce.c nlforce_d.c nlforce_p.c nlforce_s.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include "md.h"

REAL linint (REAL * y, REAL r, REAL invdr)
{

    REAL d0, d1, dm;
    REAL f0, g0, g1, g2, h1, h2, i2;
    int ic;

    d0 = r * invdr;
    ic = (int) d0;
    ic = (ic > 0) ? ic : 1;

    /* cubic interpolation using forward differences */
    d0 -= (REAL) (ic);
    d1 = (d0 - 1.0) * 0.5;
    dm = (d0 - 2.0) / 3.0;

    f0 = y[ic];
    g0 = y[ic] - y[ic - 1];
    g1 = y[ic + 1] - y[ic];
    g2 = y[ic + 2] - y[ic + 1];
    h1 = g1 - g0;
    h2 = g2 - g1;
    i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + dm * i2));

}                               /* end linint */

/******/
