

#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10
#define  EPSR   1.e-6

void gcxcpw91_spin(rmg_double_t rho_up, rmg_double_t rho_dw,
  rmg_double_t grad_up, rmg_double_t grad_dw, rmg_double_t grad, rmg_double_t *enxc,
  rmg_double_t *vxc1_up, rmg_double_t *vxc1_dw, rmg_double_t *vxc2_upup, rmg_double_t *vxc2_dwdw,
  rmg_double_t *vxc2_updw, rmg_double_t *vxc2_dwup)
{
	rmg_double_t rhotot, lim, zeta;	
        rmg_double_t sxup, sxdw, sx, v1xup, v1xdw, v2xup, v2xdw, sc, v1cup, v1cdw, v2c;


	/* exchange correction for spin up */
	if ( rho_up > SMALL && grad_up > SMALL)
		ggax( 2.0 * rho_up, 2.0 * grad_up, &sxup, &v1xup, &v2xup);
	else
	{
		sxup = 0.0;
		v1xup = 0.0;
		v2xup = 0.0;
	}

	/* exchange correction for spin dw */
	if ( rho_dw > SMALL && grad_dw > SMALL)
		ggax( 2.0 * rho_dw, 2.0 * grad_dw, &sxdw, &v1xdw, &v2xdw);
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
		zeta = (rho_up - rho_dw) / rhotot;
		if (fabs(zeta) > 1.0)
		{
			sc = 0.0;
			v1cup = 0.0;
			v1cdw = 0.0;
			v2c = 0.0;
		}
		else
		{
			lim = min(fabs(zeta), (1.0 - EPSR) );
			if (zeta >= 0.0 )
				zeta = fabs(lim);
			else if (zeta < 0.0)
				zeta = - fabs(lim);	

			if (rhotot <= SMALL || grad <= SMALL)
			{
				sc = 0.0;
				v1cup = 0.0;
				v1cdw = 0.0;
				v2c = 0.0;
			}
			else
				ggac_spin(rhotot, zeta, grad, &sc, &v1cup, &v1cdw, &v2c);
		}
		
			
	}
	else
	{
		sc = 0.0;
		v1cup = 0.0;
		v1cdw = 0.0;
		v2c = 0.0;
		
	}

	*enxc = sx + sc;
	*vxc1_up = v1xup + v1cup;
	*vxc1_dw = v1xdw + v1cdw;
	*vxc2_upup = 2.0 * v2xup + v2c;
	*vxc2_dwdw = 2.0 * v2xdw + v2c;
	*vxc2_updw = v2c;
	*vxc2_dwup = v2c;
}

