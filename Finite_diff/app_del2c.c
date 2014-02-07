/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_del2c.c *****
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
 *    rmg_double_t app_del2c(rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
 *                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
 *    Applies a 7-point finite difference discretization of the 
 *    Laplacian operator to  a matrix 
 * INPUTS
 *    a[(dimx+2)*(dimy+2)*(dimz+2)]: the matrix the Laplacian is applied to
 *    dimx, dimy, dimz: array dimentions in x, y, z direction
 *    gridhx, gridhy, gridhz:  grid spacings in crystal coordinates
 *
 *    Note: The Laplacian is not applied to the boundary layer but the 
 *          boundary layer must contain the correct values before entrance 
 *          into the function.
 * OUTPUT
 *    b[(dimx+2)*(dimy+2)*(dimz+2)]: the results of Laplacian(a)
 * PARENTS
 *    solv_pois.c, eval_residual.c
 * CHILDREN
 *    nothing
 * SOURCE
 */


#include "const.h"
#include "common_prototypes.h"
#include "fixed_dims.h"
#include "rmg_alloc.h"


#include <float.h>
#include <math.h>


rmg_double_t app_del2c (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_del2c_rmg_double(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}                               /* end app_del2c */


rmg_double_t app_del2c_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    return FD_app_del2c_rmg_float(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);

}                               /* end app_del2c */

/******/
