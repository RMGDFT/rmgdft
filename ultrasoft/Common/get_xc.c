
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_xc.c *****
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
 *   void get_xc (REAL *nrho, REAL *nrho_oppo, REAL *vxc, REAL *exc, int xctype. int spinflag)
 *   Top-level driver routine that calls a subfunction to generate
 *   a specific exchange-corrleation potential.
 * INPUTS
 *   nrho: Electronic charge density in spin-pairwised calculation, while in spin-polarized calculation,
 *        it's processor's own spin charge density (nrho already has nonlinear core charges added if any)
 *   nrho_oppo: Electronic charge density of the opposite spin in spin-polarized calculation
 *   xctype: exchange-correlation flag
 *   spinflag: flag to indicate whther do spin-polarized calculation or not
 * OUTPUT
 *   vxc: The generated exchange-corrleation potential for processor's own spin 
 *   exc: The total energy density 
 * NOTE
 *   nrho and nrho_oppo already have nonlinear core charges added if there is any in pseudopotential files
 * PARENTS
 *   get_te.c get_vxc.c
 * CHILDREN
 *   xcgga.c xcgga_spin.c xclda.c xclsda_spin.c xcgga_libxc.c xcgga_spin_libxc.c xclda_libxc.c xclsda_spin_libxc.c 	
 * SOURCE
 */


#include "main.h"
#include <float.h>
#include <math.h>



void get_xc (REAL * nrho, REAL * nrho_oppo,  REAL * vxc, REAL * exc, int xctype)
{   
	int libflag;
 
	/* libflag to indicate whether use libxc or not to calculate xc potential and energy: 
	 * libflag = 1 means using libxc; while libflag = 0 means not using libxc */

	libflag = 1;

   	/* Evaluate the XC potential and energy*/
	if (ct.spin_flag && libflag)
	{
                /* XC calculation for spin polarized case*/
    		if (xctype == LDA_PZ81)
    	       		/* Perdew Wang Physical Review B, Volume 45, Number 23, 1992 */
 			xclsda_spin_libxc (nrho, nrho_oppo, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga_spin_libxc (nrho, nrho_oppo, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");

	}
	else if ( (!ct.spin_flag) && libflag)
	{

                /* XC calculation for spin unpolarized case*/
    		if (xctype == LDA_PZ81)
 			xclda_libxc (nrho, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga_libxc (nrho, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");
	}
	else if (ct.spin_flag && (!libflag))
	{
                /* XC calculation for spin polarized case*/
    		if (xctype == LDA_PZ81)
 			xclsda_spin (nrho, nrho_oppo, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga_spin (nrho, nrho_oppo, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");

	}
	else if ( (!ct.spin_flag) && (!libflag))
	{

                /* XC calculation for spin unpolarized case*/
    		if (xctype == LDA_PZ81)
 			xclda (nrho, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga (nrho, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");
	}
}
