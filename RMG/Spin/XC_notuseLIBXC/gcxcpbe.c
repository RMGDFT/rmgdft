
#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10

void gcxcpbe (rmg_double_t rho, rmg_double_t grad, rmg_double_t * enxc, rmg_double_t * vxc1, rmg_double_t * vxc2)
{
	rmg_double_t sx, sc, v1x, v2x, v1c, v2c, arho;
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

void pbec (rmg_double_t rho, rmg_double_t grad, int iflag, rmg_double_t * sc, rmg_double_t * v1c, rmg_double_t * v2c)
{
	/* PBE correction withou LDA part
	 * iflag=0: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996) 
	 * iflag=1: J.P.Perdew et al., PRL 100, 136406 (2008)*/

	rmg_double_t ga=0.031091, be[]={0.066725, 0.046};
	rmg_double_t third, pi34, xkf, xks;
	third = 1.0 / 3.0;
	pi34 = 0.6203504908994;              /* (3 / (4 * pi))^(1 / 3) */
	xkf = 1.919158292677513;             /* (9 * pi / 4)^(1 / 3) */
	xks = 1.128379167095513;             /* sqrt (4 / pi) */

	rmg_double_t kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy;
	rmg_double_t s1, h0, dh0, ddh0;

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

void pw (rmg_double_t rs, int iflag, rmg_double_t * ec, rmg_double_t * vc)
{
	/*iflag=0: J.P.Perdew and Y. Wang, PRB 45, 13244 (1992)
	 *iflag=1: G.Ortiz and P. Ballone, PRB 50, 1391 (1994)   */
	
	rmg_double_t a=0.031091, b1=7.5957, b2=3.5876, c0=a, c1=0.046644, c2=0.00664, c3=0.01043, d0=0.4335, d1=1.4408;
	rmg_double_t lnrs, rs12, rs32, rs2, om, dom, olog;
	rmg_double_t a1[]={0.21370, 0.026481}, b3[]={1.6382, -0.46647}, b4[]={0.49294, 0.13354};

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








































































































