/************************** SVN Revision Information **************************
 **    $Id: solv_pois.c 1705 2012-04-20 12:11:31Z ebriggs $    **
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
 *   void solv_pois( rmg_float_t *vmat, rmg_float_t *fmat, rmg_float_t *work, 
 *                   int dimx, int dimy, int dimz,
 *                   rmg_float_t gridhx, rmg_float_t gridhy, rmg_float_t gridhz, rmg_float_t step)
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
 *   step: iteration step
 *   k: converts equation from poisson type to helmholtz
 * OUTPUT
 *   vmat is updated
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   app_del2c.c 
 * SOURCE
 */


#include "main.h"

void solv_pois_f (rmg_float_t * vmat, rmg_float_t * fmat, rmg_float_t * work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k)
{
    int size, idx;
    rmg_double_t scale;
    rmg_double_t diag;

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        work[idx] = ZERO;
    diag = -app_del2c_f (vmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    scale = step / diag;
    
    // Non-zero k effectively means we are solving the Helmholtz rather than Poissons equation
    if(k != 0.0) {

        for (idx = 0; idx < size; idx++)
        {

            vmat[idx] += scale * (work[idx] - k*vmat[idx] - fmat[idx]);

        }                           /* end for */

     }
     else {

        for (idx = 0; idx < size; idx++)
        {

            vmat[idx] += scale * (work[idx] - fmat[idx]);

        }                           /* end for */

     }

}                               /* end solv_pois */

/******/
