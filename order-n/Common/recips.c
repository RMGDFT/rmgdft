/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/recips.c *****
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
 *   void recips(void)
 *   generates the reciprocal lattice vectors b0, b1, b2 given the real
 *   space vectors a0, a1, a2.  the b's are units of 2pi/a.
 * INPUTS
 *   a0,a1,a2 are from ct.a0, ...
 * OUTPUT
 *   b0, b1, b2 reciprocal lattice vectors (see recips.h)
 * PARENTS
 *   init.c
 * CHILDREN
 *   cross_product.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"
#include "recips.h"

void recips(void)
{

    double mag;

    cross_product(ct.a1, ct.a2, b0);
    mag = ct.a0[0] * b0[0] + ct.a0[1] * b0[1] + ct.a0[2] * b0[2];
    b0[0] /= mag;
    b0[1] /= mag;
    b0[2] /= mag;
    cross_product(ct.a2, ct.a0, b1);
    b1[0] /= mag;
    b1[1] /= mag;
    b1[2] /= mag;
    cross_product(ct.a0, ct.a1, b2);
    b2[0] /= mag;
    b2[1] /= mag;
    b2[2] /= mag;


}                               /* end recips */

/******/
