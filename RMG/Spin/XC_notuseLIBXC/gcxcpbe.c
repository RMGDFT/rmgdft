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

void gcxcpbe (double rho, double grad, double * enxc, double * vxc1, double * vxc2)
{
	double sx, sc, v1x, v2x, v1c, v2c, arho;
	int iflag;

	iflag = 0;

	sx = v1x = v2x = 0.0;
	sc = v1c = v2c = 0.0;
	
	arho = fabs(rho);
	
	/* exchange correction */
	if (arho > SMALL )
		pbex(arho, grad, iflag, &sx, &v1x, &v2x); 

	/* correlation correction */
	if (arho > SMALL )
		pbec(arho, grad, iflag, &sc, &v1c, &v2c);

	*vxc1 = v1x + v1c;
	*vxc2 = v2x + v2c;
	*enxc = sx + sc;
}

void pbec (double rho, double grad, int iflag, double * sc, double * v1c, double * v2c)
{
	/* PBE correction withou LDA part
	 * iflag=0: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996) 
	 * iflag=1: J.P.Perdew et al., PRL 100, 136406 (2008)*/

	double ga=0.031091, be[]={0.066725, 0.046};
	double third, pi34, xkf, xks;
	third = 1.0 / 3.0;
	pi34 = 0.6203504908994;              /* (3 / (4 * pi))^(1 / 3) */
	xkf = 1.919158292677513;             /* (9 * pi / 4)^(1 / 3) */
	xks = 1.128379167095513;             /* sqrt (4 / pi) */

	double kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy;
	double s1, h0, dh0, ddh0;

	rs = pi34 / pow (rho, third);        /* rs = (3 / (4 * pi * rho))^(1 / 3) */
	pw (rs, iflag, &ec, &vc );             
	kf = xkf / rs;                      /* kf = (3 pi^2 rho)^(1/3)=(9*pi/4)^(1/3) * (1/rs) */
	ks =  xks * sqrt ( kf );
	t = grad / (2.0 * ks * rho);
	expe = exp (- ec / ga);
	af = be[iflag] / ga * ( 1.0 / (expe - 1.0));
	bf = expe * (vc - ec);
	y = af * t * t;
	xy = (1.0 + y) / (1.0 + y + y * y);
	qy = y * y * (2.0 + y) / pow ( (1.0 + y + y * y ), 2);
	s1 = 1.0 + be[iflag] / ga * t * t * xy;
	h0 = ga * log (s1);
	dh0 = be[iflag] * t * t / s1 * (- 7.0 / 3.0 * xy - qy * (af * bf / be[iflag] - 7.0 / 3.0 ) );
	ddh0 = be[iflag] / (2.0 * ks * ks * rho ) * (xy - qy ) / s1;
	*sc = h0;
	*v1c = h0 + dh0;
	*v2c = ddh0;
} 

void pw (double rs, int iflag, double * ec, double * vc)
{
	/*iflag=0: J.P.Perdew and Y. Wang, PRB 45, 13244 (1992)
	 *iflag=1: G.Ortiz and P. Ballone, PRB 50, 1391 (1994)   */
	
	double a=0.031091, b1=7.5957, b2=3.5876, c0=a, c1=0.046644, c2=0.00664, c3=0.01043, d0=0.4335, d1=1.4408;
	double lnrs, rs12, rs32, rs2, om, dom, olog;
	double a1[]={0.21370, 0.026481}, b3[]={1.6382, -0.46647}, b4[]={0.49294, 0.13354};

	/* high- and low-density formula implemented but not used in PW case
	 * (reason): inconsistencies in PBE/PW91 functionals */
	
	if ( rs < 1.0 && iflag == 1 )
	{
		/* high density formula */
		lnrs =  log (rs);
		*ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs;
		*vc = c0 * lnrs - (c1 + c0 / 3.0) + 2.0 / 3.0 * c2 * rs * lnrs - (2.0 * c3 + c2) / 3.0 * rs;
	}
	else if (rs > 100.0 && iflag == 1)
	{
		/* low density formula */
		*ec = - d0 / rs + d1 / pow (rs, 1.50 );
		*vc = - 4.0 / 3.0 * d0 / rs + 1.50 * d1 / pow (rs, 1.50);
	}
	else
	{
		/* interpolation formular */
		rs12 = sqrt (rs);
		rs32 = rs * rs12;
		rs2 = rs * rs;
		om = 2.0 * a * (b1 * rs12 + b2 * rs + b3[iflag] * rs32 + b4[iflag] * rs2);
		dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.50 * b3[iflag] * rs32 + 2.0 * b4[iflag] * rs2 );
		olog = log (1.0 + 1.0 / om );
		*ec = - 2.0 * a * (1.0 + a1[iflag] * rs ) * olog;
		*vc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1[iflag] * rs) * olog - 2.0 / 3.0 * a * (1.0 + a1[iflag] * rs ) * dom / (om * (om + 1.0 ) );
	}         /* end if */
	
}








































































































