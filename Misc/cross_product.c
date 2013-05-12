/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/cross_product.c *****
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
 *   void cross_product(rmg_double_t *a, REAL *b, REAL *c)
 *   cross product of two vectors.
 * INPUTS
 *   a: a vector with dimension of 3
 *   b: a vector with dimension of 3
 * OUTPUT
 *   c: = a x b
 * PARENTS
 *   latgen.c recips.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void cross_product (rmg_double_t * a, REAL * b, REAL * c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = -(a[0] * b[2] - a[2] * b[0]);
    c[2] = a[0] * b[1] - a[1] * b[0];
}                               /* end cross_product */

/******/
