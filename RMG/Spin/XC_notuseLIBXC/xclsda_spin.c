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
#include "common_prototypes.h"
#include <float.h>
#include <math.h>

/* crs = (9*pi/4)^(1/3) */
#define    crs    1.91915829267751281
#define    SMALL  1.e-10


void xclsda_spin(double * rho, double * rho_oppo, double * vxc, double * exc)
{
    int idx;

    double pisq3, ex, ec, vx, vx_oppo, vc, vc_oppo;
    double zeta, rs, kf;
    
    double rhox, arhox; /* the total charge in each grid point and the abosolute value of the charge */
   
    pisq3 = THREE * PI * PI;

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
	rhox = rho[idx] + rho_oppo[idx];    
        arhox = fabs(rhox);

	if (arhox < SMALL )
	{   
            vxc[idx]=0.0;
	    exc[idx]=0.0;
            continue;
        }
        
	kf = pow (pisq3 * arhox, 0.333333333333333);
        rs = crs / kf;

        zeta = (rho[idx] - rho_oppo[idx]) / arhox; 
	if ( zeta > 1.0)
		zeta = 1.0;
	else if (zeta < -1.0)
		zeta = -1.0;

	/* calculate the exchange potential and energy */
	slater_spin(arhox, zeta, &ex, &vx, &vx_oppo);

	/* calculate the correlation potential and energy */
	pz_spin(rs, zeta, &ec, &vc, &vc_oppo);
		
        /* get exchange-correlation potential and energy */
        vxc[idx] = vc + vx;
        exc[idx] = ec + ex;

    }
  	    /* end for */
}                               /* end xclsda_spin */



void slater_spin(double arhox, double zeta, double * ex, double * vx, double * vx_oppo)
{
	double vxup, vxdw, f, alpha, third, p43;
	double exup, exdw, rho13, rhoup, rhodw;

	f = -1.10783814957303361;  /* -9/8*(3/pi)^(1/3) */
	alpha = 2.0 / 3.0;
	third = 1.0 / 3.0; 
	p43 = 4.0 / 3.0;	

	rhoup = (1 + zeta) * arhox;
	rho13 = pow(rhoup, third);
	exup = f * alpha * rho13;
	vxup = p43 * f * alpha * rho13;

	rhodw = (1 - zeta) * arhox;
	rho13 = pow(rhodw, third);
	exdw = f * alpha * rho13;
	vxdw = p43 * f * alpha * rho13;

	*ex = 0.5 * ((1.0 + zeta) * exup + (1 - zeta) * exdw);
	*vx = vxup; 
	*vx_oppo = vxdw;	
	
}   /* end slater_spin */



void pz_spin(double rs, double zeta, double * ec, double * vc, double * vc_oppo)
{
	double ecu, vcu, ecp, vcp, fz, dfz;
	double p43, third;
	int iflag = 0;
	double p1, p2, p3;

	p43 = 4.0 / 3.0;
	third = 1.0 / 3.0;

	/* unpolarized part (Perdew-Zunger formula) */
	pz(rs, iflag, &ecu, &vcu);
	
	/* polarized contribution */
	pz_polarized(rs, &ecp, &vcp);

	p1 = pow((1.0 + zeta), p43);
	p2 = pow((1.0 - zeta), p43);
	p3 = pow(2.0, p43);
	fz = (p1 + p2 -2.0) / (p3 - 2.0);
	
	p1 = pow((1.0 + zeta), third);
	p2 = pow((1.0 - zeta), third);
	dfz = p43 * (p1 - p2) / (p3 - 2.0);

	*ec = ecu + fz * (ecp - ecu);
	*vc = vcu + fz * (vcp -vcu) + (ecp - ecu) * dfz * (1.0 - zeta);
	*vc_oppo = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (-1.0 -zeta);

} /* end pz_spin*/




/* J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) 
 * Spin-polarized energy and potential*/

void pz_polarized(double rs, double * ec, double * vc)
{
	double a = 0.01555, b = -0.0269, c = 0.0007, d = -0.0048;
        double gc = -0.0843, b1 = 1.3981, b2 = 0.2611;
	double lnrs, rs12, ox, dox;

	if (rs < 1.0)   /* high density formula */
	{
		lnrs = log (rs);
		*ec = a * lnrs + b + c * rs * lnrs + d * rs;
		*vc = a * lnrs + (b - a / 3.0) + 2.0 / 3.0 * c * rs * lnrs + (2.0 * d - c) / 3.0 * rs; 
	}
	else        /* interpolation formula */
	{
		rs12 = sqrt(rs);
		ox = 1.0 + b1 * rs12 + b2 * rs;
		dox = 1.0 + 7.0 / 6.0 * b1 * rs12 + 4.0 / 3.0 * b2 * rs;
		*ec = gc / ox;
		*vc = (*ec) * dox / ox;
	}      
	    /* end if */
	
}  /* end pz_polarized */





























































