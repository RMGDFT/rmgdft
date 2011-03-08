/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/solv_pois.c *****
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
 *   void solv_pois( REAL *vmat, REAL *fmat, REAL *work, 
 *                   int dimx, int dimy, int dimz,
 *                   REAL gridhx, REAL gridhy, REAL gridhz)
 *   Routines for solving poisson's equation iteratively.
 *		del**2 vmat= fmat
 *   Uses weighted Jacobi iterative method for one iteration
 * INPUTS
 *   vmat[(dimx+2)*(dimy+2)*(dimz+2)]: images must be present already
 *   fmat[(dimx+2)*(dimy+2)*(dimz+2)]: images are useless
 *   work: work space, size at least as big as vmat
 *   dimx, dimy, dimz: size of array except for images
 *   gridhx:  grid spacing in x plane
 *   gridhy:  grid spacing in y plane
 *   gridhz:  grid spacing in z plane
 * OUTPUT
 *   vmat is updated
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   app_del2c.c 
 * SOURCE
 */


#include "md.h"

void solv_pois (REAL * vmat, REAL * fmat, REAL * work,
                int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz)
{
    int size, idx;
    REAL scale;
    REAL diag;

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);

    for (idx = 0; idx < size; idx++)
        work[idx] = 0.0;
    diag = -app_del2c (vmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    scale = 1.0 / diag;
    for (idx = 0; idx < size; idx++)
    {

        vmat[idx] += scale * (work[idx] - fmat[idx]);

    }                           /* end for */

}                               /* end solv_pois */

/******/
