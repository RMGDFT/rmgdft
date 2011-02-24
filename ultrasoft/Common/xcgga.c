/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/xcgga.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
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


#include "main.h"
#include <float.h>
#include <math.h>


#define    crs    1.91915829267751281
#define    SMALL  1.e-10


void xcgga (REAL * rho, REAL * vxc, REAL * exc, int mode)
{
    int ix, iy, iz, idx, iflag;
    FP0_GRID *gx, *gy, *gz, *vgx, *vgy, *vgz, *agg, *d2rho;
    REAL d, grad, vxc1, vxc2[FP0_BASIS], enxc;
    REAL kf, pisq3, ex, vx, ec, vc, rs;


    pisq3 = THREE * PI * PI;


    /* Grab some memory */
    my_malloc (gx, 1, FP0_GRID);
    my_malloc (gy, 1, FP0_GRID);
    my_malloc (gz, 1, FP0_GRID);
    my_malloc (vgx, 1, FP0_GRID);
    my_malloc (vgy, 1, FP0_GRID);
    my_malloc (vgz, 1, FP0_GRID);
    my_malloc (agg, 1, FP0_GRID);
    my_malloc (d2rho, 1, FP0_GRID);






    /* Generate the gradient of the density */
    app_gradf (rho, gx, gy, gz);


    /* Get the Laplacian of the density */
    app6_del2f (rho, d2rho);


    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        agg->s2[idx] = sqrt (gx->s2[idx] * gx->s2[idx] +
                             gy->s2[idx] * gy->s2[idx] + gz->s2[idx] * gz->s2[idx]);

    }                           /* end for */


    /* The LDA part of exchange correlation potential and energy */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    d = fabs (rho[idx]);
	    if (d < SMALL && ct.scf_steps < 10)
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
	grad = agg->s2[idx]; 

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
    app_gradf (vxc2, vgx, vgy, vgz);


     /* add the second term gradient correction to xc potential */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	     vxc[idx] -= ( vgx->s2[idx] * gx->s2[idx] + 
	     		   vgy->s2[idx] * gy->s2[idx] + vgz->s2[idx] * gz->s2[idx] ) ;
	     vxc[idx] -= vxc2[idx] * d2rho->s2[idx];
    }





/*
    printf ("\nprint out 50_50 vxc\n");
    print_density_z_direction(50,50,vxc,FPX0_GRID,FPY0_GRID,FPZ0_GRID,20);
*/

    /* Release our memory */
    my_free (d2rho);
    my_free (agg);
    my_free (vgz);
    my_free (vgy);
    my_free (vgx);
    my_free (gz);
    my_free (gy);
    my_free (gx);

}                               /* end xcgga */
