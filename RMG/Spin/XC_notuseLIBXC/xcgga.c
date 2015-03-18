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


/****f* QMD-MGDFT/xcgga.c *****
 *
 * FUNCTION
 *   void xcgga(P0_GRID *rho, P0_GRID *vxc, P0_GRID *exc, int mode)
 *   Functions to generate the exchange-correlation potential and
 *   energy density using the generalized gradient approximation.
 * INPUTS
 *   rho: total valence charge density
 *   mode: type of GGA 
 *         0 = LDA-PZ81
 *         1 = GGA-BLYP
 *         2 = GGA-XB-CP
 *         3 = GGA-XP-CP
 *         4 = GGA-PBE
 * OUTPUT
 *   vxc: exchange correlation potential
 *   exc: exchange correlation energy  
 * PARENTS
 *   get_te.c get_vxc.c
 * CHILDREN
 *   pack_ptos.c trade_images.c app_grad.c app6_del2.c    
 * SOURCE
 */

#include "portability.h"
#include <float.h>
#include <math.h>
#include "common_prototypes.h"
#include "main.h"


#define    crs    1.91915829267751281
#define    SMALL  1.e-10


void xcgga (double * rho, double * vxc, double * exc, int mode)
{
    int idx, iflag, sizr;
    double *gx, *gy, *gz, *vgx, *vgy, *vgz, *agg, *d2rho;
    double d, grad, vxc1, *vxc2, enxc;
    double kf, pisq3, ex, vx, ec, vc, rs;
    double hxxgrid, hyygrid, hzzgrid;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    FP0_BASIS = get_FP0_BASIS();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();

    sizr = FP0_BASIS;

    pisq3 = THREE * PI * PI;


    /* Grab some memory */
    my_malloc (gx, sizr, double);
    my_malloc (gy, sizr, double);
    my_malloc (gz, sizr, double);
    my_malloc (vgx, sizr, double);
    my_malloc (vgy, sizr, double);
    my_malloc (vgz, sizr, double);
    my_malloc (agg, sizr, double);
    my_malloc (d2rho, sizr, double);


    my_malloc (vxc2, FP0_BASIS, double);




    /* Generate the gradient of the density */
    app_grad (rho, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);


    /* Get the Laplacian of the density */
    app6_del2 (rho, d2rho, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid );


    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        agg[idx] = sqrt (gx[idx] * gx[idx] +
                             gy[idx] * gy[idx] + gz[idx] * gz[idx]);

    }                           /* end for */


    /* The LDA part of exchange correlation potential and energy */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    d = fabs (rho[idx]);
	    if (d < SMALL )
	    {
		    vxc[idx] = 0.0;
		    exc[idx] = 0.0;
		    continue;
	    }

	    kf = pow (pisq3 * d, 0.333333333333333);
            rs = crs / kf;

	    /* to determine which set of monte carlo parameter to use for certain xc functionals */
	    iflag = 0;           

	    if (mode == GGA_PBE)
	    {
	            /* exchange potential and energy */
                    slater (rs, &ex, &vx); 

	            /* correlation potential and energy */
	            pw (rs, iflag , &ec, &vc);
	    }
	    else if (mode == GGA_XB_CP)
	    {
                    slater (rs, &ex, &vx); 

		    pz (rs, iflag , &ec, &vc);
	    }
	    else if (mode == GGA_XP_CP)
    {
                    slater (rs, &ex, &vx); 

	            pw (rs, iflag , &ec, &vc);
        }
	    else if (mode == GGA_BLYP)
        {
                    slater (rs, &ex, &vx); 
	
	            lyp (rs, &ec, &vc);
        }


	    vxc[idx] = vc + vx;  
	    /* local density approximation contribution to xc potential and energy */
            exc[idx] = ex + ec;
	    
    }                             /* end for */



    /* add the gradient correction for exchange correlation potential and energy */
    for (idx = 0; idx < FP0_BASIS; idx++)
        {

	d = rho[idx];
	grad = agg[idx]; 

        if (mode == GGA_PBE)
        {
	     if (fabs(d) < (1.e-6) || grad < (1.e-5)  )
	     {
		     vxc2[idx] = 0;
        }
	     else
        {
	    	     gcxcpbe (d, grad, &enxc, &vxc1, &vxc2[idx]); 

		     /* add the gradient correction to xc potential and energy now */
	             exc[idx] += enxc;
                     vxc[idx] += vxc1;    /* first term of gradient correction to potential*/
	     }	     
        }
	else if (mode == GGA_XB_CP)
        {
	     if (fabs(d) < (1.e-6) || grad < (1.e-5)  )
	     {
		     vxc2[idx] = 0;
	     }
	     else
	     {
	    	     gcxbcp (d, grad, &enxc, &vxc1, &vxc2[idx]); 
	             exc[idx] += enxc;
                     vxc[idx] += vxc1;
	     }	     
        }
        else if (mode == GGA_XP_CP)
        {
	     if (fabs(d) < (1.e-6) || grad < (1.e-5)  )
            {
		     vxc2[idx] = 0;
            }
            else
            {
	    	     gcxcpw91 (d, grad, &enxc, &vxc1, &vxc2[idx]); 
	             exc[idx] += enxc;
                     vxc[idx] += vxc1;
            }
        }
	else if (mode == GGA_BLYP)
        {
	     if (fabs(d) < (1.e-6) || grad < (1.e-5)  )
	    { 
		     vxc2[idx] = 0;
        }
        else
        {
	    	     gcxcblyp (d, grad, &enxc, &vxc1, &vxc2[idx]); 
	             exc[idx] += enxc;
                     vxc[idx] += vxc1;
	     }	     
        }

    }                           /* end for */ 



    /* Get gradient of vxc2 */
    app_grad (vxc2, vgx, vgy, vgz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);


     /* add the second term gradient correction to xc potential */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	     vxc[idx] -= ( vgx[idx] * gx[idx] + 
	     		   vgy[idx] * gy[idx] + vgz[idx] * gz[idx] ) ;
	     vxc[idx] -= vxc2[idx] * d2rho[idx];
    }





/*
    printf ("\nprint out 50_50 vxc\n");
    print_density_z_direction(50,50,vxc,FPX0_GRID,FPY0_GRID,FPZ0_GRID,20);
*/

    /* Release our memory */
    my_free (vxc2);
    my_free (d2rho);
    my_free (agg);
    my_free (vgz);
    my_free (vgy);
    my_free (vgx);
    my_free (gz);
    my_free (gy);
    my_free (gx);

}                               /* end xcgga */
