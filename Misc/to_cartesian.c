/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/to_cartesian.c *****
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
 *   void to_cartesian(rmg_double_t crystal[], rmg_double_t cartesian[])
 *   get cartesian coordinates from crystal coordinates 
 * INPUTS
 *   crystal[3]: coordinates in crystal unit
 * OUTPUT
 *   cartesian[3]: coordinates in atomic unit
 * PARENTS
 *   cdfastrlx.c fastrlx.c rmg_fastrelax.c moldyn.c get_nlop_p.c get_nlop_d.c 
 *   get_phase.c getpoi_bc.c iiforce.c init_pos.c init_wflcao.c lforce.c
 *   metric.c nlccforce.c nlforce_d.c nlforce_p.c nlforce_s.c read_data.c
 * CHILDREN
 *   nothing
 * SEE ALSO
 *   to_crystal.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "recips.h"

void to_cartesian (rmg_double_t crystal[], rmg_double_t cartesian[])
{
    int ir;

    /*  position is in crystal coordinates, get cartesian coordinates */

    for (ir = 0; ir < 3; ir++)
    {
        cartesian[ir] = crystal[0] * ct.a0[ir] + crystal[1] * ct.a1[ir] + crystal[2] * ct.a2[ir];
    }                           /* end for ir */

}                               /* end to_cartesian */

/******/
