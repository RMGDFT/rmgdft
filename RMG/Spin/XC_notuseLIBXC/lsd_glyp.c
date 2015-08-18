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



#include "portability.h"
#include <float.h>
#include <math.h>
#include "main.h"


void lsd_glyp ( double rhoa, double rhob, double grhoaa, double grhoab2, double grhobb, double * sc, double * v1ca, double * v2ca, double * v1cb, double * v2cb, double * v2cab )
{
	double a=0.04918, b=0.132, c=0.2533, d=0.349;
	double rho, rm3, dr, or, dor, der, dder;
	double dlaa, dlab, dlbb, dlaaa, dlaab, dlaba, dlabb, dlbba, dlbbb;

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



void lsd_lyp ( double rho, double zeta, double * elyp, double * valyp, double * vblyp )
{
	/* C. Lee, W. Yang and R. G. Parr, PRB 37, 785 ( 1988 ) 
	 * only the correlation for LDA part */

	double rhoa, rhob, rm3, dr, en1, or, dor, en2, de1a, de1b, de2a, de2b;
	double small_val=1.0e-24, a=0.04918, b=0.132, c=0.2533, d=0.349, cf=2.87123400018819108;

	rhoa = rho * 0.5 * ( 1.0 + zeta );
	rhoa = rmg_max ( rhoa, small_val );
	rhob = rho * 0.5 * ( 1.0 - zeta );
	rhob = rmg_max ( rhob, small_val );

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












































































