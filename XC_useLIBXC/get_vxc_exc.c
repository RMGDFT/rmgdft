#include "main.h"
#include <float.h>
#include <math.h>



void get_vxc_exc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype)
{   
	int libflag;
 
	/* libflag to indicate whether use libxc or not to calculate xc potential and energy: 
	 * libflag = 1 means using libxc; while libflag = 0 means not using libxc */

	libflag = 0;

   	/* Evaluate the XC potential and energy*/
	if (ct.spin_flag && libflag)
	{
                /* XC calculation for spin polarized case*/
    		if (xctype == LDA_PZ81 || (xctype == MGGA_TB09 ))
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
    		if (xctype == LDA_PZ81 || (xctype == MGGA_TB09))
 			xclda_libxc (nrho, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga_libxc (nrho, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");
	}
	else if (ct.spin_flag && (!libflag))
	{
                /* XC calculation for spin polarized case*/
    		if (xctype == LDA_PZ81  || (xctype == MGGA_TB09))
 			xclsda_spin (nrho, nrho_oppo, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga_spin (nrho, nrho_oppo, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");

	}
	else if ( (!ct.spin_flag) && (!libflag))
	{

                /* XC calculation for spin unpolarized case*/
    		if (xctype == LDA_PZ81  || (xctype == MGGA_TB09))
 			xclda (nrho, vxc, exc);
		else if ( (xctype == GGA_BLYP) || (xctype == GGA_XB_CP) || (xctype == GGA_XP_CP) || (xctype == GGA_PBE) )
			 xcgga (nrho, vxc, exc, xctype);
		else
        		error_handler ("Unknown exchange-correlation functional");
	}
}
