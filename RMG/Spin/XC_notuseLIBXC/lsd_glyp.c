

#include "main.h"
#include <float.h>
#include <math.h>


void lsd_glyp ( rmg_double_t rhoa, rmg_double_t rhob, rmg_double_t grhoaa, rmg_double_t grhoab2, rmg_double_t grhobb, rmg_double_t * sc, rmg_double_t * v1ca, rmg_double_t * v2ca, rmg_double_t * v1cb, rmg_double_t * v2cb, rmg_double_t * v2cab )
{
	rmg_double_t a=0.04918, b=0.132, c=0.2533, d=0.349;
	rmg_double_t rho, rm3, dr, or, dor, der, dder;
	rmg_double_t dlaa, dlab, dlbb, dlaaa, dlaab, dlaba, dlabb, dlbba, dlbbb;

	rho = rhoa + rhob;
	rm3 = pow( rho, - 1.0 / 3.0 );
	dr = 1.0 + d * rm3;
	or = exp ( - c * rm3 ) / dr * pow ( rm3, 11 );
	dor = - 1.0 / 3.0 * pow ( rm3, 4 ) * or * ( 11.0 / rm3 - c - d / dr );
	der = c * rm3 + d * rm3 / dr;
	dder = 1.0 / 3.0 * ( d * d * pow ( rm3, 5 ) / dr / dr - der / rho );
	dlaa = - a * b * or * ( rhoa * rhob / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhoa / rho ) - rhob * rhob ); 
	dlab = - a * b * or * ( rhoa * rhob / 9.0 * ( 47.0 - 7.0 * der ) - 4.0 / 3.0 * rho * rho );
	dlbb = - a * b * or * ( rhoa * rhob / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhob / rho ) - rhoa * rhoa );

	dlaaa = dor / or * dlaa - a * b * or * ( rhob / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhoa / rho ) - rhoa * rhob / 9.0 * ( (  3.0 + rhoa / rho) * dder + ( der - 11.0 ) * rhob / rho / rho ) );
	
	dlaab = dor / or * dlaa - a * b * or * ( rhoa / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhoa / rho ) - rhoa * rhob / 9.0 * ( ( 3.0 + rhoa / rho ) * dder - ( der - 11.0 ) * rhoa / rho / rho ) - 2.0 * rhob );
	
	dlaba = dor / or * dlab - a * b * or * ( rhob / 9.0 * ( 47.0 - 7.0 * der ) - 7.0 / 9.0 * rhoa * rhob * dder - 8.0 / 3.0 * rho );
	
	dlabb = dor / or * dlab - a * b * or * ( rhoa / 9.0 * ( 47.0 - 7.0 * der ) - 7.0 / 9.0 * rhoa * rhob * dder - 8.0 / 3.0 * rho );
	
	dlbba = dor / or * dlbb - a * b * or * ( rhob / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhob / rho ) - rhoa * rhob / 9.0 * ( ( 3.0 + rhob / rho ) * dder - ( der - 11.0 ) * rhob / rho / rho ) - 2.0 * rhoa );
	
	dlbbb = dor / or * dlbb - a * b * or * ( rhoa / 9.0 * ( 1.0 - 3 * der - ( der - 11.0 ) * rhob /rho ) - rhoa * rhob / 9.0 * ( ( 3.0 + rhob / rho ) * dder + ( der - 11.0 ) * rhoa / rho / rho ) );

	*sc = dlaa * grhoaa * grhoaa + dlab * grhoab2 + dlbb * grhobb * grhobb;
	*v1ca = dlaaa * grhoaa * grhoaa + dlaba * grhoab2 + dlbba * grhobb * grhobb;
	*v1cb = dlaab * grhoaa * grhoaa + dlabb * grhoab2 + dlbbb * grhobb * grhobb;
	*v2ca = 2.0 * dlaa;
	*v2cb = 2.0 * dlbb;
	*v2cab = dlab;

	*sc = (*sc) / rho;

}



void lsd_lyp ( rmg_double_t rho, rmg_double_t zeta, rmg_double_t * elyp, rmg_double_t * valyp, rmg_double_t * vblyp )
{
	/* C. Lee, W. Yang and R. G. Parr, PRB 37, 785 ( 1988 ) 
	 * only the correlation for LDA part */

	rmg_double_t rhoa, rhob, rm3, dr, en1, or, dor, en2, de1a, de1b, de2a, de2b;
	rmg_double_t small=1.e-24, a=0.04918, b=0.132, c=0.2533, d=0.349, cf=2.87123400018819108;

	rhoa = rho * 0.5 * ( 1.0 + zeta );
	rhoa = max ( rhoa, small );
	rhob = rho * 0.5 * ( 1.0 - zeta );
	rhob = max ( rhob, small );

	rm3 = pow ( rho, - 1.0 / 3.0 );
	dr = 1.0 + d * rm3;
	en1 = 4.0 * a * rhoa * rhob / rho / dr;
	or = exp ( -c * rm3 ) / dr * pow ( rm3, 11 );
	dor = - 1.0 / 3.0 * pow ( rm3, 4 ) * or * ( 11.0 / rm3 - c - d / dr );
	en2 = pow ( 2.0, 11.0 / 3.0 ) * cf * a * b * or * rhoa * rhob * ( pow ( rhoa, 8.0 / 3.0 ) + pow ( rhob, 8.0 / 3.0 ) );
	*elyp = - ( en1 + en2 ) / rho;
	de1a = - en1 * ( 1.0 / 3.0 * d * pow ( rm3, 4 ) / dr + 1.0 / rhoa - 1.0 / rho );
	de1b = - en1 * ( 1.0 / 3.0 * d * pow ( rm3, 4 ) / dr + 1.0 / rhob - 1.0 / rho );

	de2a = - pow ( 2.0, 11.0 / 3.0 ) * cf * a * b * ( dor * rhoa * rhob * ( pow ( rhoa, 8.0 / 3.0 ) + pow ( rhob, 8.0 / 3.0 ) ) 
			+ or * rhob * ( 11.0 / 3.0 * pow ( rhoa, 8.0 / 3.0 ) + pow ( rhob, 8.0 / 3.0 ) ) );
	
	de2b = - pow ( 2.0, 11.0 / 3.0 ) * cf * a * b * ( dor * rhoa * rhob * ( pow ( rhoa, 8.0 / 3.0 ) + pow ( rhob, 8.0 / 3.0 ) ) 
			+ or * rhoa * ( 11.0 / 3.0 * pow ( rhob, 8.0 / 3.0 ) + pow ( rhoa, 8.0 / 3.0 ) ) );
	*valyp = de1a + de2a;
	*vblyp = de1b + de2b;
	
}












































































