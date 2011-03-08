/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/get_mehr.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *  void get_mehr(void)
 *    Generates the mehrstellen coefficients.  
 * INPUTS
 *   nothing
 * OUTPUT
 *  
 * PARENTS
 *   init.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "md.h"

void get_mehr(void)
{

    REAL ihx, ihy, ihz;

    ct.Bc = 0.5;
    ct.Bx = ct.By = ct.Bz = 1.0 / 12.0;

    ihx = 1. / (ct.hxgrid * ct.hxgrid * ct.xside * ct.xside);
    ihy = 1. / (ct.hygrid * ct.hygrid * ct.yside * ct.yside);
    ihz = 1. / (ct.hzgrid * ct.hzgrid * ct.zside * ct.zside);

    ct.Ac = (-4. / 3.) * (ihx + ihy + ihz);

    ct.Ax = (5.0 / 6.0) * ihx + (ct.Ac / 8.0);
    ct.Ay = (5.0 / 6.0) * ihy + (ct.Ac / 8.0);
    ct.Az = (5.0 / 6.0) * ihz + (ct.Ac / 8.0);

    ct.Axy = (1.0 / 12.0) * (ihx + ihy);
    ct.Axz = (1.0 / 12.0) * (ihx + ihz);
    ct.Ayz = (1.0 / 12.0) * (ihy + ihz);

    return;

}                               /* end get_mehr */

/********/
