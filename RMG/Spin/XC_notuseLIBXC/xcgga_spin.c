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

/* crs = (9*pi/4)^(1/3) */
#define    crs    1.91915829267751281
#define    SMALL  1.e-15
#define     ETA   1.e-10

/* up and dn are convenient for naming the processor's own spin the the opposite spin, 
 doesn't mean the real physics */

void xcgga_spin(double * rho_up, double * rho_dw, double * vxc_up, double * exc, int mode)

{
    int idx, sizr;
    double *d2rho_up, *d2rho_dw, *agg, *agg_updw2;
    double *gx_up, *gy_up, *gz_up, *agg_up;
    double *gx_dw, *gy_dw, *gz_dw, *agg_dw; 
    double *gx_vuu, *gy_vuu, *gz_vuu, *gx_vud, *gy_vud, *gz_vud; 

    double pisq3, ex, ec, vxup, vxdw, vcup, vcdw;
    double rhotot, arhox, zeta, rs, kf;
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

    double *vxc2_upup, *vxc2_updw, vxc2_dwdw, vxc2_dwup;
    double vxc1_up, vxc1_dw, grad_up, grad_dw, grad, grad_updw2, enxc, gx, gy, gz;


    /* Grab some memory */ 
    /* to hold gradient of spin up charge density */ 
    my_malloc (gx_up, sizr, double);
    my_malloc (gy_up, sizr, double);
    my_malloc (gz_up, sizr, double);

    /* to hold gradient of spin down charge density */
    my_malloc (gx_dw, sizr, double);
    my_malloc (gy_dw, sizr, double);
    my_malloc (gz_dw, sizr, double);

    /* to hold the absolute of the gradient of total, up and down density */
    my_malloc (agg, sizr, double);
    my_malloc (agg_up, sizr, double);
    my_malloc (agg_dw, sizr, double);

    if ( mode == GGA_BLYP )
    	my_malloc (agg_updw2, sizr, double);
        /* to holde  (grad rhoup) \dot (grad rhodw)  */

    
    /* to hold laplaciant of the spin up and down charge density */
    my_malloc (d2rho_up,  sizr, double);
    my_malloc (d2rho_dw,  sizr, double);

    /* to hold the gradient of potentials */
    my_malloc (gx_vuu,  sizr, double);
    my_malloc (gy_vuu,  sizr, double);
    my_malloc (gz_vuu,  sizr, double);

    my_malloc (gx_vud,  sizr, double);
    my_malloc (gy_vud,  sizr, double);
    my_malloc (gz_vud,  sizr, double);
    my_malloc (vxc2_upup, sizr, double);
    my_malloc (vxc2_updw, sizr, double);


    /* Generate the gradient of the density */
    app_grad (rho_up, gx_up, gy_up, gz_up, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);    /* spin up density */
    app_grad (rho_dw, gx_dw, gy_dw, gz_dw, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);    /* spin down density */
    

    /* Get the Laplacian of the density */
    app6_del2 (rho_up, d2rho_up, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid );
    app6_del2 (rho_dw, d2rho_dw, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid );

    
    
    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	gx = gx_up[idx] + gx_dw[idx];
	gy = gy_up[idx] + gy_dw[idx];
	gz = gz_up[idx] + gz_dw[idx]; 

        agg[idx] = sqrt ( gx * gx + gy * gy + gz * gz);
        
	/*Absolute value of grad(rho_up)*/    
        agg_up[idx] = sqrt (gx_up[idx] * gx_up[idx] +
                             gy_up[idx] * gy_up[idx] + gz_up[idx] * gz_up[idx]);
        
	/*Absolute value of grad(rho_dw) */
        agg_dw[idx] = sqrt (gx_dw[idx] * gx_dw[idx] +
                             gy_dw[idx] * gy_dw[idx] + gz_dw[idx] * gz_dw[idx]);
        
	/*  (grad rhoup) \dot (grad rhodw)  */
	if (mode == GGA_BLYP)
        	agg_updw2[idx] = gx_up[idx] * gx_dw[idx] +
                	             gy_up[idx] * gy_dw[idx] + gz_up[idx] * gz_dw[idx];

    }                           
   


    /* Caculate the LSDA part of exchange correlation potential and energy */

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	rhotot = rho_up[idx] + rho_dw[idx];    
        arhox = fabs(rhotot);

	if (arhox < SMALL )
	{   
            vxc_up[idx]=0.0;
	    exc[idx]=0.0;
            continue;
        }
        
	kf = pow (pisq3 * arhox, 0.333333333333333);
        rs = crs / kf;

        zeta = (rho_up[idx] - rho_dw[idx]) / arhox; 
	if ( zeta > 1.0)
		zeta = 1.0;
	else if (zeta < -1.0)
		zeta = -1.0;

	if (mode == GGA_PBE)
	{
		/* calculate the exchange potential and energy */
		slater_spin(arhox, zeta, &ex, &vxup, &vxdw);
	
		/* calculate the correlation potential and energy */
		pw_spin(rs, zeta, &ec, &vcup, &vcdw);
	}
	else if (mode == GGA_XB_CP)
	{
		slater_spin(arhox, zeta, &ex, &vxup, &vxdw);
		
		pz_spin(rs, zeta, &ec, &vcup, &vcdw);

	}
	else if (mode == GGA_XP_CP)
	{
		slater_spin(arhox, zeta, &ex, &vxup, &vxdw);

		pw_spin(rs, zeta, &ec, &vcup, &vcdw);
	}
	else if (mode == GGA_BLYP)
	{
		slater_spin(arhox, zeta, &ex, &vxup, &vxdw);

		lsd_lyp (arhox, zeta, &ec, &vcup, &vcdw);
	}

		
        vxc_up[idx] = vcup + vxup;
        exc[idx] = ec + ex;

    }
  	    


    
    /* Add the gradient correction for exchange correlation potential and energy */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

	grad_up = agg_up[idx];
	grad_dw = agg_dw[idx];
	grad = agg[idx];
	    

        if (mode == GGA_PBE)
        {
		gcxcpbe_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
        }                      
	else if (mode == GGA_XB_CP)
        {
		
		gcxbcp_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
		
        }                      
	else if (mode == GGA_XP_CP)
        {
		
		gcxcpw91_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
		
        }                       
	else if (mode == GGA_BLYP)
        {
		
		grad_updw2 = agg_updw2[idx];

		gcxcblyp_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad_updw2, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
		
        }                       

    }                           /* end for */ 


		
    /* get the exchange-correlation potential correction
     * vxcup += vxc1up - div( vxc2upup * grad(rhoup) + vxc2updw * grad(rhodw) ) 
     * vxcdw += vxc1dw - div( vxc2dwdw * grad(rhodw) + vxc2dwup * grad(rhoup) )*/ 

    /* Generate the gradient of the second term exchange-correlation potential vxc2*/
    app_grad (vxc2_upup, gx_vuu, gy_vuu, gz_vuu, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);
    app_grad (vxc2_updw, gx_vud, gy_vud, gz_vud, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);




    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    vxc_up[idx] -= ( gx_vuu[idx] * gx_up[idx] + 
		           gy_vuu[idx] * gy_up[idx] + gz_vuu[idx] * gz_up[idx] );
	    vxc_up[idx] -=  vxc2_upup[idx] * d2rho_up[idx] ;
	    
	    vxc_up[idx] -= ( gx_vud[idx] * gx_dw[idx] + 
		           gy_vud[idx] * gy_dw[idx] + gz_vud[idx] * gz_dw[idx] );
	    vxc_up[idx] -= vxc2_updw[idx] * d2rho_dw[idx];
    }


    /* Release our memory */  
    my_free (vxc2_updw);
    my_free (vxc2_upup);
    my_free (gx_up);
    my_free (gy_up);
    my_free (gz_up);

    my_free (gx_dw);
    my_free (gy_dw);
    my_free (gz_dw);

    my_free (agg);
    my_free (agg_up);
    my_free (agg_dw); 

    my_free(d2rho_up);
    my_free(d2rho_dw);

    my_free (gx_vuu);
    my_free (gy_vuu);
    my_free (gz_vuu);

    my_free (gx_vud);
    my_free (gy_vud);
    my_free (gz_vud);

    if (mode == GGA_BLYP)
    	my_free (agg_updw2);
    

}                               /* end xcgga */
