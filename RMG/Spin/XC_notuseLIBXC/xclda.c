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


#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>

#define    crs    1.91915829267751281
#define    SMALL  1.e-10


void xclda(double * rho, double * vxc, double * exc)
{
    int  idx, iflag;
    double d, kf;
    
    double pisq3, ex, vx, ec, vc;
    double rs;
   
    pisq3 = THREE * PI * PI;

    /* Now get the potential */
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
	d = fabs(rho[idx]);    /* use the absolute value of rho since intepolated charge may be negative*/    
        
	if (d < SMALL )
	{   
            vxc[idx]=0.0;
	    exc[idx]=0.0;
            continue;
        }
        
	kf = pow (pisq3 * d, 0.333333333333333);
        rs = crs / kf;
	
	/* exchange potential and energy */
        slater(rs, &ex, &vx); 
	
	/* correlation potential and energy */
	iflag = 0; 
	pz(rs, iflag , &ec, &vc);
        
	vxc[idx] = vc + vx;
        exc[idx] = ex + ec;

    }

}                               /* end xclda */ 


/* slater exchange with alpha = 2.0 / 3.0 */

void slater(double rs, double * ex, double * vx)
{
	double f, alpha;
	f = -0.687247939924714;   /* -9/8*(3/(2*pi))^(2/3) */
        alpha = 2.0 / 3.0;
	*ex = f * alpha / rs;
	*vx = 4.0 / 3.0 * f * alpha / rs;	
}



/* LDA parameterization from Monte Carlo data 
 * iflag=0: J.P.Perdew and A. Zunger, PRB 23, 5048 (1981)
 * iflag=1: G.Ortiz and P.Ballone, PRB 50, 1391 (1994)*/

void pz(double rs, int iflag, double * ec, double * vc)
{
	double lnrs, rs12, ox, dox;
	double a[]={0.0311, 0.031091}, b[]={-0.048, -0.046644};
	double c[]={0.0020, 0.00419}, d[]={-0.0116, -0.00983};
	double gc[]={-0.1423, -0.103756}, b1[]={1.0529, 0.56371}, b2[]={0.3334, 0.27358};

	if ( rs < 1.0)         /* high density formula*/
	{
		lnrs = log(rs);
		*ec = a[iflag] * lnrs + b[iflag] + c[iflag] * rs * lnrs + d[iflag] * rs;
		*vc = a[iflag] * lnrs + (b[iflag] - a[iflag] / 3.0) + 2.0 / 3.0 * c[iflag] * rs * lnrs + (2.0 * d[iflag] -c[iflag]) / 3.0 * rs;

	}
	else        /* interpolation formula */
	{
		rs12 = sqrt(rs);
		ox = 1.0 + b1[iflag] * rs12 + b2[iflag] * rs;
		dox = 1.0 + 7.0 / 6.0 * b1[iflag] * rs12 + 4.0 / 3.0 * b2[iflag] * rs;
		*ec = gc[iflag] / ox;
		*vc = (*ec) * dox / ox;
	}
}
	
























