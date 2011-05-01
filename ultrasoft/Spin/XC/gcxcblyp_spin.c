

#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10
#define  EPSR   1.e-6

void gcxcblyp_spin (REAL rho_up, REAL rho_dw,
  REAL grad_up, REAL grad_dw, REAL grad_updw2, REAL *enxc,
  REAL *vxc1_up, REAL *vxc1_dw, REAL *vxc2_upup, REAL *vxc2_dwdw,
  REAL *vxc2_updw, REAL *vxc2_dwup)
{
	REAL rhotot;	
        REAL sxup, sxdw, sx, v1xup, v1xdw, v2xup, v2xdw;
	REAL sc, v1cup, v1cdw, v2cup, v2cdw, v2cupdw;

	
	/* exchange correction for spin up */
	if ( rho_up > SMALL && grad_up > SMALL)
		becke88_spin ( rho_up,  grad_up, &sxup, &v1xup, &v2xup);
	else
	{
		sxup = 0.0;
		v1xup = 0.0;
		v2xup = 0.0;
	}

	/* exchange correction for spin dw */
	if ( rho_dw > SMALL && grad_dw > SMALL)
		becke88_spin ( rho_dw, grad_dw, &sxdw, &v1xdw, &v2xdw);
	else
	{
		sxdw = 0.0;
		v1xdw = 0.0;
		v2xdw = 0.0;
	}


        rhotot = rho_up + rho_dw; 
	
	/* average correction of spin up and down */
	if (rhotot > SMALL)
		sx =  (sxup * rho_up  + sxdw * rho_dw) / rhotot;
	else
		sx =  (sxup * rho_up  + sxdw * rho_dw) / SMALL;


	/*correlation correction */
	if (rhotot > EPSR )
	{
		if (rhotot < 1.e-20)
		{
			sc = 0.0;
			v1cup = 0.0;
			v1cdw = 0.0;
			v2cup = 0.0;
			v2cdw = 0.0;
			v2cupdw = 0.0;
		}
		else
		{
			lsd_glyp (rho_up, rho_dw, grad_up, grad_updw2, grad_dw, &sc, &v1cup, &v1cdw, &v2cup, &v2cdw, &v2cupdw);
		}
		
			
	}
	else
	{
		sc = 0.0;
		v1cup = 0.0;
		v1cdw = 0.0;
		v2cup = 0.0;
		v2cdw = 0.0;
		v2cupdw = 0.0;
		
	}

	*enxc = sx + sc;
	*vxc1_up = v1xup + v1cup;
	*vxc1_dw = v1xdw + v1cdw;
	*vxc2_upup = v2xup + v2cup;
	*vxc2_dwdw = v2xdw + v2cdw;
	*vxc2_updw = v2cupdw;
	*vxc2_dwup = v2cupdw;
	
}

