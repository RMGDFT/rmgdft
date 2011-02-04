/************************** SVN Revision Information **************************
 **    $Id: get_vxc_soft.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/get_vxc_soft.c *****
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
 *   void get_vxc(REAL *rho, REAL *rhocore, REAL *vxc)
 *   Top-level driver routine that calls a subfunction to generate
 *   a specific exchange-corrleation potential.
 * INPUTS
 *   rho: Electronic charge density
 *   rhocore: Core charge density if Non-linear core corrections are being used.
 *   vxc: The generated exchange-corrleation potential.
 * OUTPUT
 *   vxc: The generated exchange-corrleation potential.
 * PARENTS
 *   init.c scf.c
 * CHILDREN
 *   xclda_pz81.c xcgga.c xclsd_pz81.c
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>



void get_vxc_soft (REAL * rho, REAL * rhocore, REAL * vxc)
{

    int idx, ispin, item;
    REAL *nrho, esum, *exc;
    int nn, ione = 1;

    nn = ct.vh_pbasis;          /*shuchun wang */
    esum = ddot (&nn, vxc, &ione, rho, &ione);
    ct.Evxcold_rho = ct.vel_f * real_sum_all (esum);


    /*begin shuchun wang */
    if (ct.spin)
        my_malloc_init( nrho, ct.spin * ct.vh_pbasis, REAL );
    if (!ct.spin)
        my_malloc_init( nrho, ct.vh_pbasis, REAL );
    my_malloc_init( exc, ct.vh_pbasis, REAL );

    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        item = ct.spin * ct.vh_pbasis;
        for (idx = 0; idx < ct.vh_pbasis; idx++)
            nrho[idx + item] = rho[idx + item] + rhocore[idx];
    }
    if (ct.spin)
        xclsd_pz81 (nrho, vxc);
    if (!ct.spin)
    {
        switch (ct.xctype)
        {
        case LDA_PZ81:         /* LDA Perdew Zunger 81 */

            xclda_pz81 (nrho, vxc);
            break;

        case GGA_BLYP:         /* GGA X-Becke C-Lee Yang Parr */

            xcgga (nrho, vxc, exc, ct.xctype);
            break;

        case GGA_XB_CP:        /* GGA X-Becke C-Perdew */

            xcgga (nrho, vxc, exc, ct.xctype);
            break;

        case GGA_XP_CP:        /* GGA X-Perdew C-Perdew */

            xcgga (nrho, vxc, exc, ct.xctype);
            break;

        case GGA_PBE:          /* GGA Perdew, Burke, Ernzerhof */

            xcgga (nrho, vxc, exc, ct.xctype);
            break;

        default:
            error_handler ("Unknown exchange-correlation functional");

        }                       /* end switch */
    }                           /*end if */

    my_free(nrho);
    my_free(exc);
    /* end shuchun wang */

}                               /* end get_vxc */
