/************************** SVN Revision Information **************************
 **    $Id: get_vxc.c 1066 2009-08-31 18:41:09Z froze $    **
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



void get_vxc_spin (REAL * rho_f, REAL * rho_oppo,  REAL * rhocore_f, REAL * vxc_f)
{

    int idx;
    REAL *exc, *nrho_up, *nrho_dn;    
    
     /* up and down is convenient for naming the processor's own spin and the opposite spin,  
        but not the real meaning of spin up and down */

    my_malloc (exc, 3 * FP0_BASIS, REAL);
    nrho_up = exc + FP0_BASIS;
    nrho_dn = exc + 2 * FP0_BASIS;

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        nrho_up[idx] = rho_f[idx];
	nrho_dn[idx] = rho_oppo[idx];
	
    }

    /* Add in the core charges */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        nrho_up[idx] += rhocore_f[idx] / 2.0;         /*Divide by two since the core charges are spin up and spin down paired*/
	nrho_dn[idx] += rhocore_f[idx] / 2.0;
    }



    /* Evaluate the XC potential */
    switch (ct.xctype)
    {
    case LDA_PZ81:       /* Perdew Wang Physical Review B, Volume 45, Number 23, 1992 */
        

        xclsda_spin (nrho_up, nrho_dn, vxc_f, exc);
        break;

    case GGA_BLYP:             /* GGA X-Becke C-Lee Yang Parr */

        xcgga_spin (nrho_up, nrho_dn, vxc_f, exc, ct.xctype);
        break;

    case GGA_XB_CP:            /* GGA X-Becke C-Perdew */

        xcgga_spin (nrho_up, nrho_dn, vxc_f, exc, ct.xctype);
        break;

    case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */

        xcgga_spin (nrho_up, nrho_dn, vxc_f, exc, ct.xctype);
        break;

    case GGA_PBE:              /* GGA Perdew, Burke, Ernzerhof */

        xcgga_spin (nrho_up, nrho_dn, vxc_f, exc, ct.xctype);
        break;

    default:
        error_handler ("Unknown exchange-correlation functional");

    }                           /* end switch */

    my_free (exc);

}                               /* end get_vxc */



/******/
