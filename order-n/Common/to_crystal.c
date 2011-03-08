/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/to_crystal.c *****
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
 *   void to_crystal(REAL crystal[], REAL cartesian[])
 *   get crystal coordinates from cartesian coordinates 
 * INPUTS
 *   cartesian[3]: coordinates in atomic unit
 * OUTPUT
 *   crystal[3]: coordinates in crystal unit
 * PARENTS
 *   cdfastrlx.c fastrlx.c md_fastrelax.c moldyn.c symmetry.c
 * CHILDREN
 *   nothing
 * SEE ALSO
 *   to_cartesian.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"
#include "recips.h"

#define     SR3         1.732050807569

void to_crystal(REAL crystal[], REAL cartesian[])
{
    int ir;

    if (ct.ibrav == HEXAGONAL)
    {

        crystal[0] = (cartesian[0] - cartesian[1] / SR3) / ct.celldm[0];
        crystal[1] = cartesian[1] / (SR3 / 2.0) / ct.celldm[0];
        crystal[2] = cartesian[2] * b2[2];

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else if (ct.ibrav == CUBIC_PRIMITIVE)
    {

        crystal[0] = cartesian[0] / ct.celldm[0];
        crystal[1] = cartesian[1] / ct.celldm[0];
        crystal[2] = cartesian[2] / ct.celldm[0];

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else if (ct.ibrav == ORTHORHOMBIC_PRIMITIVE)
    {

        crystal[0] = cartesian[0] / ct.celldm[0];
        crystal[1] = cartesian[1] / (ct.celldm[0] * ct.celldm[1]);
        crystal[2] = cartesian[2] / (ct.celldm[0] * ct.celldm[2]);

        if (crystal[0] < 0.0)
            crystal[0] += 1.0;
        if (crystal[1] < 0.0)
            crystal[1] += 1.0;
        if (crystal[2] < 0.0)
            crystal[2] += 1.0;
        if (crystal[0] > 1.0)
            crystal[0] -= 1.0;
        if (crystal[1] > 1.0)
            crystal[1] -= 1.0;
        if (crystal[2] > 1.0)
            crystal[2] -= 1.0;

    }
    else
    {

        for (ir = 0; ir < 3; ir++)
        {
            crystal[ir] = cartesian[0] * b0[ir] + cartesian[1] * b1[ir] + cartesian[2] * b2[ir];
        }                       /* end for ir */

    }                           /* end if */


}                               /* end to_crystal  */

/******/
