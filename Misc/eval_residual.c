/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/eval_residual.c *****
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
 *   void eval_residual( rmg_double_t *mat, REAL *f_mat, int dimx, int dimy, int dimz, 
 *                       rmg_double_t gridhx, REAL gridhy, REAL gridhz, REAL *res )
 *  Evaluate the residual of Poisson equation 
 *      f_mat - (del**2)mat 
 * INPUTS
 *   mat: array dimensioned (dimx+2)*(dimy+2)*(dimz+2)
 *	  (+2 for images of neighboring cells)
 *	  (images must be present already)
 *   f_mat: array with same dimensions as mat
 *   dimx, dimy, dimz: size of array except for images
 *   gridhx: grid spacing in x plane
 *   gridhy: grid spacing in y plane
 *   gridhz: grid spacing in z plane
 * OUTPUT
 *   res: result of f_mat - (del**2) mat
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   app_del2c.c
 * SOURCE
 */


#include "main.h"


void eval_residual (rmg_double_t * mat, REAL * f_mat, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, REAL gridhy, REAL gridhz, REAL * res)
{
    int size, idx;

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        res[idx] = 0.0;
    app_del2c (mat, res, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    for (idx = 0; idx < size; idx++)
        res[idx] = f_mat[idx] - res[idx];


}                               /* end eval_residual */

/******/
