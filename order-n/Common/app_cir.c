/************************** SVN Revision Information **************************
 **    $Id: app_cir.c 587 2006-08-18 22:57:56Z miro $    **
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
 *    REAL app_cir(REAL *a, REAL *b, int dimx, int dimy, int dimz)
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

#include "md.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>


void app_cir(REAL * a, REAL * b, int dimx, int dimy, int dimz)
{

#if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc();
#endif

    switch (ct.ibrav)
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:
        app_cir_ortho(a, b, dimx, dimy, dimz);
        break;

    case CUBIC_BC:
        app_cir_bcc(a, b, dimx, dimy, dimz);
        break;

    case CUBIC_FC:
        app_cir_fcc(a, b, dimx, dimy, dimz);
        break;

    case HEXAGONAL:
        app_cir_hex(a, b, dimx, dimy, dimz);
        break;

    default:
        error_handler("Lattice type not implemented");

    }                           /* end switch */



#if MD_TIMERS
    time2 = my_crtc();
    rmg_timings(APPCIR_TIME, (time2 - time1), 0);
#endif

}                               /* end app_cir */


/******/
