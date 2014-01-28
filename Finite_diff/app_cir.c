/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_cir.c *****
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
 *    rmg_double_t app_cir(rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz)
 *    Top level driver routine for applying the Mehrstellen RHS operator.
 * INPUTS
 *    a[(dimx+2) * (dimy+2) * (dimz+2)]: the matrix to be applied
 *    dimx, dimy, dimz: array dimentions in x, y, z direction
 * OUTPUT
 *    b[dimx * dimy * dimz]:  RHS(a)
 * PARENTS
 *    get_vh.c init_wf.c mg_eig_state.c subdiag_mpi.c subdiag_smp.c 
 * CHILDREN
 *    app_cir_ortho.c app_cir_bcc.c app_cir_fcc.c app_cir_hex.c
 * SOURCE
 */

#include "make_conf.h"
#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>


void app_cir (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz)
{


    switch (ct.ibrav)
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:
        app_cir_ortho (a, b, dimx, dimy, dimz);
        break;

    case CUBIC_BC:
        app_cir_bcc (a, b, dimx, dimy, dimz);
        break;

    case CUBIC_FC:
        app_cir_fcc (a, b, dimx, dimy, dimz);
        break;

    case HEXAGONAL:
        app_cir_hex (a, b, dimx, dimy, dimz);
        break;

    default:
        error_handler ("Lattice type not implemented");

    }                           /* end switch */



}                               /* end app_cir */


/******/
