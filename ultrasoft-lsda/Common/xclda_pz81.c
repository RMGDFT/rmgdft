/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/xclda_pz81.c *****
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
 *   void xclda_pz81(REAL *rho, REAL *vxc)
 *   void exclda_pz81(REAL *rho, REAL *exc)
 *   Functions to generate the exchange-correlation potential and
 *   energy density using the Ceperley - Alder exchange-correlation
 *   potential and energy as parameterized by Perdew and Zunger, 
 *   Phys. Rev. B23, 5048 (1981)
 * INPUTS
 *   rho: total valence charge density
 * OUTPUT
 *   vxc: exchange correlation potential
 *   exc: exchange correlation energy
 * PARENTS
 *   get_te.c get_vxc.c 
 * CHILDREN
 *   nothing 
 * SOURCE
 */

#include "main.h"
#include <float.h>
#include <math.h>

#define MA  -0.4582E0
#define MB  -0.0480E0
#define MC   0.0311E0
#define MD   0.0020E0
#define ME  -0.0116E0
#define NB  -0.1423E0
#define NC   1.0529E0
#define ND   0.3334E0



/* Computes exchange correlation energy density */
/* rho is the normalized RS electron charge density */
double e_pz (double rho)
{

    double rs, eval;

    if (rho <= 1.0e-15)
    {

        return 0.0;

    }                           /* end if */


    rs = (3.0 / (rho * 4.0 * PI));
    rs = pow (rs, 0.3333333333333333333333333e0);
    eval = MA / rs;

    if (rs < 1.0)
    {

        eval = eval + MB + (MC + MD * rs) * log (rs) + ME * rs;

    }
    else
    {

        eval = eval + NB / (1.0 + NC * sqrt (rs) + ND * rs);

    }                           /* end if */


    return eval;

}                               /* end e_pz */





/*  Compute exchange-correlation potential (in Hartrees).   */
/*  Subroutine adapted from Hamann                          */
/*  rho is the normalized real-space electron charge density */

double mu_pz (double rho)
{

    double rs, sqrs, eval;

    if (rho <= 1.0e-15)
    {

        return 0.0;

    }                           /* end if */


    rs = (3.0 / (rho * 4.0 * PI));
    rs = pow (rs, 0.3333333333333333333333333);
    eval = (4.0 * MA / 3.0) / rs;


    if (rs < 1.0)
    {

        eval = eval + (MB - MC / 3.0) + (2.0 * ME - MD) * rs / 3.0 +
            (MC + 2.0 / 3.0 * MD * rs) * log (rs);

    }
    else
    {

        sqrs = sqrt (rs);
        eval = eval + NB *
            (1.0 + 7.0 / 6.0 * NC * sqrs + 4.0 * ND * rs / 3.0) /
            ((1.0 + NC * sqrs + ND * rs) * (1.0 + NC * sqrs + ND * rs));

    }                           /* end if */

    return eval;

}                               /* end mu_pz */




void xclda_pz81 (REAL * rho, REAL * vxc_f)
{

    int i;

    for (i = 0; i < FP0_BASIS; i++)
    {

        vxc_f[i] = mu_pz (rho[i]);

    }                           /* end for */

    /*my_barrier(); */


}                               /* end xclda_pz81 */




void exclda_pz81 (REAL * rho, REAL * exc)
{

    int i;

    for (i = 0; i < FP0_BASIS; i++)
    {

        exc[i] = e_pz (rho[i]);

    }                           /* end for */

    /*my_barrier(); */


}                               /* end exclda_pz81 */

/******/
