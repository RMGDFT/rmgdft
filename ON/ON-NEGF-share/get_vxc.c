/************************* SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/get_vxc.c *****
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
 *   void get_vxc(rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *vxc)
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


#include "main.h"
#include "prototypes_on.h"
#include <float.h>
#include <math.h>



void get_vxc(rmg_double_t * rho, rmg_double_t * rhocore, rmg_double_t * vxc)
{

    int idx, ispin, item;
    rmg_double_t *nrho, esum, *exc;
    int nn, ione = 1;

    nn = ct.vh_pbasis;          /*shuchun wang */
    esum = ddot(&nn, vxc, &ione, rho, &ione);
    ct.Evxcold_rho = get_vel_f() * real_sum_all(esum, pct.grid_comm);


    /*begin shuchun wang */
    if (ct.spin)
        my_malloc_init( nrho, ct.spin * ct.vh_pbasis, rmg_double_t );
    if (!ct.spin)
        my_malloc_init( nrho, ct.vh_pbasis, rmg_double_t );
    my_malloc_init( exc, ct.vh_pbasis, rmg_double_t );

    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        item = ct.spin * ct.vh_pbasis;
        for (idx = 0; idx < ct.vh_pbasis; idx++)
            nrho[idx + item] = rho[idx + item] + rhocore[idx];
    }
    if (ct.spin)
        xclsd_pz81(nrho, vxc);
    if (!ct.spin)
    {
        switch (ct.xctype)
        {
        case LDA_PZ81:         /* LDA Perdew Zunger 81 */

            xclda_pz81(nrho, vxc, nn );
            break;

        case GGA_BLYP:         /* GGA X-Becke C-Lee Yang Parr */

            xcgga(nrho, vxc, exc, ct.xctype);
            break;

        case GGA_XB_CP:        /* GGA X-Becke C-Perdew */

            xcgga(nrho, vxc, exc, ct.xctype);
            break;

        case GGA_XP_CP:        /* GGA X-Perdew C-Perdew */

            xcgga(nrho, vxc, exc, ct.xctype);
            break;

        case GGA_PBE:          /* GGA Perdew, Burke, Ernzerhof */

            xcgga(nrho, vxc, exc, ct.xctype);
            break;

        default:
            error_handler("Unknown exchange-correlation functional");

        }                       /* end switch */
    }                           /*end if */

    my_free(nrho);
    my_free(exc);
    /* end shuchun wang */

}                               /* end get_vxc */
