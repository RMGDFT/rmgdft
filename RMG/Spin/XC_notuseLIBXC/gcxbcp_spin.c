/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10
#define  EPSR   1.e-6

void gcxbcp_spin (double rho_up, double rho_dw,
  double grad_up, double grad_dw, double grad, double *enxc,
  double *vxc1_up, double *vxc1_dw, double *vxc2_upup, double *vxc2_dwdw,
  double *vxc2_updw, double *vxc2_dwup)
{
	double rhotot, lim, zeta;	
        double sxup, sxdw, sx, v1xup, v1xdw, v2xup, v2xdw, sc, v1cup, v1cdw, v2c;


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
			lim = rmg_min(fabs(zeta), (1.0 - EPSR) );
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
				perdew86_spin(rhotot, zeta, grad, &sc, &v1cup, &v1cdw, &v2c);
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
	*vxc2_upup = v2xup + v2c;
	*vxc2_dwdw = v2xdw + v2c;
	*vxc2_updw = v2c;
	*vxc2_dwup = v2c;
}

