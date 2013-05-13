

#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10

void gcxcpw91 (rmg_double_t rho, rmg_double_t grad, rmg_double_t * enxc, rmg_double_t * vxc1, rmg_double_t * vxc2)
{
	rmg_double_t sx, sc, v1x, v2x, v1c, v2c, arho;
	int iflag;

	iflag = 0;

	sx = v1x = v2x = 0.0;
	sc = v1c = v2c = 0.0;
	
	arho = fabs(rho);
	
	/* exchange correction */
	if (arho > SMALL )
		ggax(arho, grad, &sx, &v1x, &v2x); 

	/* correlation correction */
	if (arho > SMALL )
		ggac(arho, grad, &sc, &v1c, &v2c);

	*vxc1 = v1x + v1c;
	*vxc2 = v2x + v2c;
	*enxc = sx + sc;
}

void ggax(rmg_double_t rho, rmg_double_t grho, rmg_double_t * sx, rmg_double_t * v1x, rmg_double_t * v2x)
{
	/* Perdew-Wang GGA (PW91), exchange part: J.P.Perdew et al., PRB 46, 6671 (1992) */
	
	rmg_double_t f1=0.19645, f2=7.7956, f3=0.2743, f4=0.1508, f5=0.004;
	rmg_double_t fp1=-0.019292021296426, fp2=0.161620459673995;     
	/* fp1=-3/(16 pi)*(3 pi^2)^(-1/3)   fp2=(1/2)(3 pi^2)^(-1/3)*/
	rmg_double_t rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, dbs, dls;

	rhom43 = pow ( rho, - 4.0/3.0 );
	s = fp2 * grho * rhom43;
	s2 = s * s;
	s3 = s2 * s;
	s4 = s3 * s;
	exps = f4 * exp ( -100.0 * s2 );
	as = f3 - exps - f5 * s2;
	sa2b8 = sqrt ( 1.0 + f2 * f2 * s2 );
	shm1 = log ( f2 * s + sa2b8 );
	bs = 1.0 + f1 * s * shm1 + f5 * s4;
	das = ( 200.0 * exps - 2.0 * f5 ) * s;
	dbs = f1 * ( shm1 + f2 * s / sa2b8 ) + 4.0 * f5 * s3;
	dls = das / as - dbs / bs;
	*sx = fp1 * grho * grho * rhom43 * as / bs;
	*v1x = - 4.0 / 3.0 * (*sx) / rho * ( 1.0 + s * dls );
	*v2x = fp1 * rhom43 * as / bs * ( 2.0 + s * dls );
	*sx = (*sx) / rho;
}

void ggac (rmg_double_t rho, rmg_double_t grho, rmg_double_t * sc, rmg_double_t * v1c, rmg_double_t * v2c)
{
	/* Perdew-Wang GGA (PW91), correlation part: J.P.Perdew et al., PRB 46, 6671 (1992) */
	
	rmg_double_t al=0.09, pa=0.023266, pb=7.389e-6, pc=8.723, pd=0.472;
	rmg_double_t cx=-0.001667, cxc0=0.002568, cc0=-cx + cxc0;
	rmg_double_t third=1.0/3.0, pi34=0.6203504908994; /*pi34 = (3/(4*pi) )^(1/3)*/
	rmg_double_t nu=15.755920349483144, be=nu*cc0, xkf=1.919158292677513, xks=1.128379167095513;
	/* nu=(16/pi)*(3 pi^2)^(1/3), xkf=(9 pi/4)^(1/3), xks=sqrt(4/pi) */

	rmg_double_t kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy, qy, s1;
	rmg_double_t h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, ddh1;
	int iflag;
	
	rs = pi34 / pow ( rho, third );
        rs2 = rs * rs;
	rs3 = rs2 * rs;

	iflag=0;
	pw (rs, iflag, &ec, &vc);
	kf = xkf / rs;
	ks = xks * sqrt (kf);
	t = grho / ( 2.0 * ks * rho);
	expe = exp ( - 2.0 * al * ec / ( be * be) );
	af = 2.0 * al / be * ( 1.0 / ( expe - 1.0 ) );
	bf = expe * ( vc - ec );
	y = af * t * t;
	xy = ( 1.0 + y ) / ( 1.0 + y + y * y );
	qy = y * y * ( 2.0 + y ) / pow ( ( 1.0 + y + y * y ), 2 );
	s1 = 1.0 + 2.0 * al / be * t * t * xy;
	h0 = be * be / ( 2.0 * al ) * log ( s1 );
	dh0 = be * t * t / s1 * ( - 7.0 / 3.0 * xy - qy * ( af * bf / be -7.0 / 3.0) );
	ddh0 = be / ( 2.0 * ks * ks * rho ) * ( xy - qy ) / s1;
	ee = - 100.0 * pow ( ( ks / kf * t ), 2 );
	cna = cxc0 + pa * rs + pb * rs2;
	dcna = pa * rs + 2.0 * pb * rs2;
	cnb = 1.0 + pc * rs + pd * rs2 + (1.e4) * pb * rs3;
	dcnb = pc * rs + 2.0 * pd * rs2 + (3.e4) * pb * rs3;
	cn = cna / cnb - cx;
	dcn = dcna / cnb - cna * dcnb / ( cnb * cnb );
	h1 = nu * ( cn - cc0 - 3.0 / 7.0 * cx ) * t * t * exp ( ee );
	dh1 = - third * ( h1 * ( 7.0 + 8.0 * ee ) + nu * t * t * exp ( ee ) * dcn );
	ddh1 = 2.0 * h1 * ( 1.0 + ee ) * rho / ( grho * grho );

	*sc = ( h0 + h1 );
	*v1c = h0 + h1 + dh0 + dh1;
	*v2c = ddh0 + ddh1;
		
	
}



void ggac_spin (rmg_double_t rho, rmg_double_t zeta, rmg_double_t grho, rmg_double_t * sc, rmg_double_t * v1cup, rmg_double_t * v1cdw, rmg_double_t * v2c)
{
	/* Perdew-Wang GGA (PW91), correlation part: J.P.Perdew et al., PRB 46, 6671 (1992) -Spin polarized */

	rmg_double_t al=0.09, pa=0.023266, pb=7.389e-6, pc=8.723, pd=0.472;
	rmg_double_t cx=-0.001667, cxc0=0.002568, cc0=-cx + cxc0;
	rmg_double_t third=1.0/3.0, pi34=0.6203504908994; /*pi34 = (3/(4*pi) )^(1/3)*/
	rmg_double_t nu=15.755920349483144, be=nu*cc0, xkf=1.919158292677513, xks=1.128379167095513;
	/* nu=(16/pi)*(3 pi^2)^(1/3), xkf=(9 pi/4)^(1/3), xks=sqrt(4/pi) */

	rmg_double_t kf, ks, rs, rs2, rs3, ec, vcup, vcdw, t, expe, af, y, xy, qy, s1;
	rmg_double_t h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, ddh1;
	rmg_double_t fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, dh0zdw, dh1zup, dh1zdw;
	
	rs = pi34 / pow ( rho, third );
        rs2 = rs * rs;
	rs3 = rs2 * rs;

	pw_spin (rs, zeta, &ec, &vcup, &vcdw);
	kf = xkf / rs;
	ks = xks * sqrt (kf);
	
	fz = pow ( ( 1.0 + zeta ), 2.0 / 3.0 ) + pow ( ( 1.0 - zeta ), 2.0 / 3.0 );
	fz = 0.5 * fz;
	fz2 = fz * fz;
	fz3 = fz2 * fz;
	fz4 = fz3 * fz;
	dfz = pow ( ( 1.0 + zeta ), - 1.0 / 3.0 ) - pow ( ( 1.0 - zeta ), - 1.0 / 3.0 );
	dfz = dfz / 3.0;
	
	t = grho / ( 2.0 * fz * ks * rho);
	expe = exp ( - 2.0 * al * ec / ( fz3 * be * be) );
	af = 2.0 * al / be * ( 1.0 / ( expe - 1.0 ) );
	bfup = expe * ( vcup - ec ) / fz3;
	bfdw = expe * ( vcdw - ec ) / fz3;
	y = af * t * t;
	xy = ( 1.0 + y ) / ( 1.0 + y + y * y );
	qy = y * y * ( 2.0 + y ) / pow ( ( 1.0 + y + y * y ), 2 );
	s1 = 1.0 + 2.0 * al / be * t * t * xy;
	h0 = fz3 * be * be / ( 2.0 * al ) * log ( s1 );
	dh0up = be * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * ( af * bfup / be -7.0 / 3.0) );
	dh0dw = be * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * ( af * bfdw / be -7.0 / 3.0) );
	
	dh0zup = ( 3.0 * h0 / fz - be * t * t * fz2 / s1 * ( 2.0 * xy - qy * ( 3.0 * af * expe * ec / fz3 / be + 2.0 ) ) ) * dfz * ( 1.0 - zeta );
	dh0zdw = - ( 3.0 * h0 / fz - be * t * t * fz3 / s1 * ( 2.0 * xy - qy * ( 3.0 * af * expe * ec / fz3 / be + 2.0 ) ) ) * dfz * ( 1.0 + zeta );
	
	ddh0 = be * fz / ( 2.0 * ks * ks * rho ) * ( xy - qy ) / s1;
	ee = - 100.0 * fz4 * pow ( ( ks / kf * t ), 2 );
	cna = cxc0 + pa * rs + pb * rs2;
	dcna = pa * rs + 2.0 * pb * rs2;
	cnb = 1.0 + pc * rs + pd * rs2 + (1.e4) * pb * rs3;
	dcnb = pc * rs + 2.0 * pd * rs2 + (3.e4) * pb * rs3;
	cn = cna / cnb - cx;
	dcn = dcna / cnb - cna * dcnb / ( cnb * cnb );
	h1 = nu * ( cn - cc0 - 3.0 / 7.0 * cx ) * fz3 * t * t * exp ( ee );
	dh1 = - third * ( h1 * ( 7.0 + 8.0 * ee ) + fz3 * nu * t * t * exp ( ee ) * dcn );
	ddh1 = 2.0 * h1 * ( 1.0 + ee ) * rho / ( grho * grho );

	dh1zup = ( 1.0 - zeta ) * dfz * h1 * ( 1.0 + 2.0 * ee / fz );
	dh1zdw = - ( 1.0 + zeta ) * dfz * h1 * ( 1.0 + 2.0 * ee / fz );
	
	*sc = ( h0 + h1 );
	*v1cup = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup;
	*v1cdw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw;
	*v2c = ddh0 + ddh1;
		
	
}










































