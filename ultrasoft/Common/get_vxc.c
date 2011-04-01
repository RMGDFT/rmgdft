/************************** SVN Revision Information **************************
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
 *   xclda_pz81.c xcgga.c
 * SOURCE
 */


#include "main.h"
#include <float.h>
#include <math.h>



void get_vxc (REAL * rho_f, REAL * rhocore_f, REAL * vxc_f)
{

    int idx;
    REAL *exc, *nrho;


    my_malloc (exc, 2 * FP0_BASIS, REAL);
    nrho = exc + FP0_BASIS;

    for (idx = 0; idx < FP0_BASIS; idx++)
        nrho[idx] = rho_f[idx];

    /* Add in the core charges */
    for (idx = 0; idx < FP0_BASIS; idx++)
        nrho[idx] += rhocore_f[idx];


    switch (ct.xctype)
    {
    case LDA_PZ81:             /* LDA Perdew Zunger 81 */

        /* xclda_pz81 (nrho, vxc_f); */

	/* incoporate both the Perdew Zunger 1981 and Ortiz Ballone 1994, default is PZ 1981 */
        //xclda (nrho, vxc_f, exc);
        xclda_libxc (nrho, vxc_f, exc);
        break;

    case GGA_BLYP:             /* GGA X-Becke C-Lee Yang Parr */

        //xcgga (nrho, vxc_f, exc, ct.xctype);
        xcgga_libxc (nrho, vxc_f, exc, ct.xctype);
        break;

    case GGA_XB_CP:            /* GGA X-Becke C-Perdew */

        //xcgga (nrho, vxc_f, exc, ct.xctype);
        xcgga_libxc (nrho, vxc_f, exc, ct.xctype);
        break;

    case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */

        //xcgga (nrho, vxc_f, exc, ct.xctype);
        xcgga_libxc (nrho, vxc_f, exc, ct.xctype);
        break;

    case GGA_PBE:              /* GGA Perdew, Burke, Ernzerhof */

        //xcgga (nrho, vxc_f, exc, ct.xctype);
        xcgga_libxc (nrho, vxc_f, exc, ct.xctype);
        break;

    default:
        error_handler ("Unknown exchange-correlation functional");

    }                           /* end switch */


    my_free (exc);


}                               /* end get_vxc */



/******/
