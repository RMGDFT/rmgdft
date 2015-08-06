/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


/****f* QMD-MGDFT/xclda_pz81.c *****
 *
 * FUNCTION
 *   void xclda_pz81(double *rho, double *vxc)
 *   void exclda_pz81(double *rho, double *exc)
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

#include <float.h>
#include <math.h>
#include "const.h"
#include "rmg_xc.h"

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




void xclda_pz81 (double * rho, double * vxc_f, int n)
{

    int i;

    for (i = 0; i < n; i++)
    {

        vxc_f[i] = mu_pz (rho[i]);

    }                           /* end for */

    /*my_barrier(); */


}                               /* end xclda_pz81 */




void exclda_pz81 (double * rho, double * exc, int n)
{

    int i;

    for (i = 0; i < n; i++)
    {

        exc[i] = e_pz (rho[i]);

    }                           /* end for */

    /*my_barrier(); */


}                               /* end exclda_pz81 */

/******/
