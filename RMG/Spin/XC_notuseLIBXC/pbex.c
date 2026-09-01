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


void pbex (double rho, double grho, int iflag, double * sx, double * v1x, double * v2x)
{
	/* PBE exchange without Slater exchange */
	/* iflag=0 J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
	 * iflag=1 "revised" PBE: Y.Zhang et al., PRL 80, 890 (1998)
	 * iflag=2 PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008) */
	
	double kf, s1, s2, ds, dsg, exunif, fx;
	double dxunif, dfx, f1, f2, f3, dfx1;
	double third, c1, c2, c5;
	third = 1.0 / 3.0;
	c1 = 0.75 / PI;
	c2 = 3.093667726280136; /* (3 pi^2)^(1/3) */
	c5 = 4.0 / 3.0;
	
	double k[]={0.804, 1.2450, 0.804}, mu[]={0.21951, 0.21951, 0.12345679012345679012};

	kf = c2 * pow (rho, third);    /* kf = (3 pi^2 |rho|)^(1/3)*/
	dsg = 0.5 / kf;  
        s1 = grho * dsg / rho;     /* s = |grho| / (2 * kf * |rho|) */	
	s2 = s1 * s1;
	ds = -c5 * s1;   

	/* correction of exchane energy */
	f1 = s2 * mu[iflag] / k[iflag];
	f2 = 1.0 + f1;
	f3 = k[iflag] / f2;
	fx = k[iflag] - f3;     /* fx = k - k / (1 + us^2 / k)*/
	exunif = -c1 * kf;     /* - (3 / 4) * (3 / pi) ^(1 / 3) * rho^(1 / 3) */
	*sx = exunif * fx;     /* exchange energy correction */

	/* correction of exchange potential */
	dxunif = exunif * third;
	dfx1 =  f2 * f2;
	dfx = 2.0 * mu[iflag] * s1 / dfx1;
	*v1x = (*sx) + dxunif * fx + exunif * dfx * ds;
	*v2x = exunif * dfx * dsg / grho;
}


void pbec_spin (double rho, double zeta, double grho, int iflag, double * sc, double * v1cup, double * v1cdw, double * v2c)
{
	/* PBE correction without LDA part - spin-polarized */
	/* iflag=1 J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996) 
	 * iflag=2 J.P.Perdew et al., PRL 100, 136406 (2008)*/

	double ga=0.031091, be[]={0.066725, 0.046};
	double third, third2, pi34, xkf, xks;
	third = 1.0 / 3.0;
	third2 = 2.0 / 3.0;
	pi34 = 0.6203504908994;    /* (3/(4 *pi))^(1/3) */
	xkf = 1.919158292677513;    /* (9*pi/4)^(1/3) */
	xks = 1.128379167095513;	/* sqrt(pi/4)*/

	double kf, ks, rs, ec, vcup, vcdw, t, expe, af, y, xy, qy, s1, h0, ddh0;
	double fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, dh0zdw; 

	rs = pi34 / pow (rho, third); 
	pw_spin(rs, zeta, &ec, &vcup, &vcdw);
	kf = xkf / rs;
	ks = xks * sqrt (kf);
	fz = pow ( (1.0 + zeta), third2 );
	fz += pow ( (1.0 - zeta), third2 );
	fz = 0.5 * fz;
	fz2 = fz * fz;
	fz3 = fz2 * fz;
	fz4 = fz3 * fz;
	dfz = pow ( (1.0 + zeta), (-third) );
	dfz -= pow ( (1.0 - zeta), (-third) );
	dfz = dfz / 3.0;
	t = grho / (2.0 * fz * ks * rho);
	expe = exp ( - ec / (fz3 * ga) );
	af = be[iflag] / ga * ( 1.0 / (expe - 1.0) );	
	bfup = expe * (vcup - ec) / fz3;
	bfdw = expe * (vcdw - ec) / fz3;
	y = af * t * t;
	xy = (1.0 + y) / (1.0 + y + y * y);
	qy = y * y * (2.0 + y ) / pow ( (1.0 + y + y * y), 2);
        s1 = 1.0 + be[iflag] / ga * t * t * xy;
	h0 = fz3 * ga * log (s1);
	dh0up = be[iflag] * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * (af * bfup / be[iflag] - 7.0 / 3.0) );
	dh0dw = be[iflag] * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * (af * bfdw / be[iflag] - 7.0 / 3.0) );
	dh0zup = ( 3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * ( 2.0 * xy - qy * ( 3.0 * af * expe * ec / fz3 / be[iflag] + 2.0 ) ) ) * dfz * (1.0 - zeta);	
	dh0zdw = - (3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * ( 2.0 * xy - qy * ( 3.0 * af * expe * ec / fz3 / be[iflag] + 2.0 ) )  ) * dfz * (1.0 + zeta);
	ddh0 = be[iflag] * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1;
	*sc = h0;
	*v1cup = h0 + dh0up + dh0zup;
	*v1cdw = h0 + dh0dw + dh0zdw;
	*v2c = ddh0;
}



void pw_spin (double rs, double zeta, double * ec, double * vcup, double * vcdw)
{
	/* J.P.Perdew and Y. Wang, PRB 45, 13244 (1992) */

	/* xc parameter, unpolarized */
	 double a=0.031091, a1=0.21370, b1=7.5957, b2=3.5876, b3=1.6382, b4=0.49294;

	/* xc parameter, polarized */
	double ap=0.015545, a1p=0.20548, b1p=14.1189, b2p=6.1977, b3p=3.3662, b4p=0.62517;

	/* xc parameter, antiferro */
	double aa=0.016887, a1a=0.11125, b1a=10.357, b2a=3.6231, b3a=0.88026, b4a=0.49671;

	double fz0=1.709921;

	double rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz;
	double om, dom, olog, epwc, vpwc;
	double omp, domp, ologp, epwcp, vpwcp;
	double oma, doma, ologa, alpha, vpwca;
	double third, third4;
	third = 1.0 / 3.0;
	third4 = 4.0 / 3.0;

	/* if (rs < 0.5) then high density formular
	 * else if (rs > 100.0) then low density formular
	 * else interpolation formular 
	 * (not implemented yet) */

	zeta2 = zeta * zeta;
	zeta3 = zeta2 * zeta;
	zeta4 = zeta3 * zeta;
	rs12 = sqrt(rs);
	rs32 = rs * rs12;
	rs2 = rs * rs;

	/* unpolarized */
	om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2);
	dom = 2.0 * a * ( 0.5 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 + 2.0 * b4 * rs2);
	olog = log ( 1.0 + 1.0 / om );
	epwc = - 2.0 * a * ( 1.0 + a1 * rs ) * olog;
	vpwc = - 2.0 * a * ( 1.0 + 2.0 / 3.0 * a1 * rs ) * olog - 2.0 / 3.0 * a * ( 1.0 + a1 * rs ) * dom / ( om * ( om + 1.0 ) );

	/* polarized */
	omp = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2);
	domp = 2.0 * ap * ( 0.5 * b1p * rs12 + b2p * rs + 1.5 * b3p * rs32 + 2.0 * b4p * rs2);
	ologp = log ( 1.0 + 1.0 / omp );
	epwcp = - 2.0 * ap * ( 1.0 + a1p * rs ) * ologp;
	vpwcp = - 2.0 * ap * ( 1.0 + 2.0 / 3.0 * a1p * rs ) * ologp - 2.0 / 3.0 * ap * ( 1.0 + a1p * rs ) * domp / ( omp * ( omp + 1.0 ) );
     			     
	 /* antiferro */     
	oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2);
	doma = 2.0 * aa * ( 0.5 * b1a * rs12 + b2a * rs + 1.5 * b3a * rs32 + 2.0 * b4a * rs2);
	ologa = log ( 1.0 + 1.0 / oma );
	alpha =  2.0 * aa * ( 1.0 + a1a * rs ) * ologa;
	vpwca = + 2.0 * aa * ( 1.0 + 2.0 / 3.0 * a1a * rs ) * ologa + 2.0 / 3.0 * aa * ( 1.0 + a1a * rs ) * doma / ( oma * ( oma + 1.0 ) ); 

	fz = ( pow ( (1.0 + zeta), third4 ) + pow ( (1.0 - zeta), third4 ) - 2.0 ) / ( pow (2.0, third4) - 2.0 ); 
	dfz = ( pow ( (1.0 + zeta), third ) - pow ( (1.0 - zeta), third ) ) * 4.0 / ( 3.0 * ( pow ( 2.0, third4 ) - 2.0 ) ); 

	*ec = epwc + alpha * fz * (1.0 - zeta4) / fz0 + (epwcp - epwc) * fz * zeta4;
	*vcup = vpwc + vpwca * fz * (1.0 -zeta4) / fz0 + (vpwcp - vpwc) * fz * zeta4 + ( alpha / fz0 * ( dfz * (1.0 - zeta4) - 4.0 * fz * zeta3 ) + (epwcp - epwc ) * ( dfz * zeta4 + 4.0 * fz * zeta3 ) ) * (1.0 - zeta);

	*vcdw = vpwc + vpwca * fz * (1.0 -zeta4) / fz0 + (vpwcp - vpwc) * fz * zeta4 - ( alpha / fz0 * ( dfz * (1.0 - zeta4) - 4.0 * fz * zeta3 ) + (epwcp - epwc ) * ( dfz * zeta4 + 4.0 * fz * zeta3 ) ) * (1.0 + zeta);

}


































































