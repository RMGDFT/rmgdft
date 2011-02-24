/************************** SVN Revision Information **************************
 **    $Id: xcgga.c 1142 2010-09-13 21:30:11Z btan $    **
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

/* crs = (9*pi/4)^(1/3) */
#define    crs    1.91915829267751281
#define    SMALL  1.e-15
#define     ETA   1.e-10

/* up and dn are convenient for naming the processor's own spin the the opposite spin, 
 doesn't mean the real physics */

void xcgga_spin(REAL * rho_up, REAL * rho_dw, REAL * vxc_up, REAL * exc, int mode)

{
#if 1
    int idx;
    FP0_GRID *d2rho_up, *d2rho_dw, *agg, *agg_updw2;
    FP0_GRID *gx_up, *gy_up, *gz_up, *agg_up;
    FP0_GRID *gx_dw, *gy_dw, *gz_dw, *agg_dw; 
    FP0_GRID *gx_vuu, *gy_vuu, *gz_vuu, *gx_vud, *gy_vud, *gz_vud; 

    REAL pisq3, ex, ec, vxup, vxdw, vcup, vcdw;
    REAL rhotot, arhox, zeta, rs, kf;
    
    pisq3 = THREE * PI * PI;

    REAL vxc2_upup[FP0_BASIS], vxc2_updw[FP0_BASIS], vxc2_dwdw, vxc2_dwup;
    REAL vxc1_up, vxc1_dw, grad_up, grad_dw, grad, grad_updw2, enxc, gx, gy, gz;


    /* Grab some memory */ 
    /* to hold gradient of spin up charge density */ 
    my_malloc (gx_up, 1, FP0_GRID);
    my_malloc (gy_up, 1, FP0_GRID);
    my_malloc (gz_up, 1, FP0_GRID);

    /* to hold gradient of spin down charge density */
    my_malloc (gx_dw, 1, FP0_GRID);
    my_malloc (gy_dw, 1, FP0_GRID);
    my_malloc (gz_dw, 1, FP0_GRID);

    /* to hold the absolute of the gradient of total, up and down density */
    my_malloc (agg, 1, FP0_GRID); 
    my_malloc (agg_up, 1, FP0_GRID);
    my_malloc (agg_dw, 1, FP0_GRID);

    if ( mode == GGA_BLYP )
    	my_malloc (agg_updw2, 1, FP0_GRID); 
        /* to holde  (grad rhoup) \dot (grad rhodw)  */

    
    /* to hold laplaciant of the spin up and down charge density */
    my_malloc (d2rho_up, 1, FP0_GRID);
    my_malloc (d2rho_dw, 1, FP0_GRID); 

    /* to hold the gradient of potentials */
    my_malloc (gx_vuu, 1, FP0_GRID);
    my_malloc (gy_vuu, 1, FP0_GRID);
    my_malloc (gz_vuu, 1, FP0_GRID);

    my_malloc (gx_vud, 1, FP0_GRID);
    my_malloc (gy_vud, 1, FP0_GRID);
    my_malloc (gz_vud, 1, FP0_GRID);


    /* Generate the gradient of the density */
    app_gradf (rho_up, gx_up, gy_up, gz_up);    /* spin up density */
    app_gradf (rho_dw, gx_dw, gy_dw, gz_dw);    /* spin down density */
    

    /* Get the Laplacian of the density */
    app6_del2f (rho_up, d2rho_up);
    app6_del2f (rho_dw, d2rho_dw);

    
    
    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	gx = gx_up->s2[idx] + gx_dw->s2[idx];
	gy = gy_up->s2[idx] + gy_dw->s2[idx];
	gz = gz_up->s2[idx] + gz_dw->s2[idx]; 

        agg->s2[idx] = sqrt ( gx * gx + gy * gy + gz * gz);
        
	/*Absolute value of grad(rho_up)*/    
        agg_up->s2[idx] = sqrt (gx_up->s2[idx] * gx_up->s2[idx] +
                             gy_up->s2[idx] * gy_up->s2[idx] + gz_up->s2[idx] * gz_up->s2[idx]);
        
	/*Absolute value of grad(rho_dw) */
        agg_dw->s2[idx] = sqrt (gx_dw->s2[idx] * gx_dw->s2[idx] +
                             gy_dw->s2[idx] * gy_dw->s2[idx] + gz_dw->s2[idx] * gz_dw->s2[idx]);
        
	/*  (grad rhoup) \dot (grad rhodw)  */
	if (mode == GGA_BLYP)
        	agg_updw2->s2[idx] = gx_up->s2[idx] * gx_dw->s2[idx] +
                	             gy_up->s2[idx] * gy_dw->s2[idx] + gz_up->s2[idx] * gz_dw->s2[idx];

    }                           
   


    /* Caculate the LSDA part of exchange correlation potential and energy */

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	rhotot = rho_up[idx] + rho_dw[idx];    
        arhox = fabs(rhotot);

	if (arhox < SMALL && ct.scf_steps < 10)
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

	grad_up = agg_up->s2[idx];
	grad_dw = agg_dw->s2[idx];
	grad = agg->s2[idx];
	    

        if (mode == GGA_PBE)
        {
#if 1		
		gcxcpbe_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
#endif		
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
		
		grad_updw2 = agg_updw2->s2[idx];

		gcxcblyp_spin (rho_up[idx], rho_dw[idx], grad_up, grad_dw, grad_updw2, &enxc, 
				&vxc1_up, &vxc1_dw, &vxc2_upup[idx], &vxc2_dwdw, &vxc2_updw[idx], &vxc2_dwup);

		exc[idx] += enxc;
		vxc_up[idx] += vxc1_up;
		
        }                       

    }                           /* end for */ 


		
    /* get the exchange-correlation potential correction
     * vxcup += vxc1up - div( vxc2upup * grad(rhoup) + vxc2updw * grad(rhodw) ) 
     * vxcdw += vxc1dw - div( vxc2dwdw * grad(rhodw) + vxc2dwup * grad(rhoup) )*/ 

#if 1
    /* Generate the gradient of the second term exchange-correlation potential vxc2*/
    app_gradf (vxc2_upup, gx_vuu, gy_vuu, gz_vuu);    
    app_gradf (vxc2_updw, gx_vud, gy_vud, gz_vud);   




    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    vxc_up[idx] -= ( gx_vuu->s2[idx] * gx_up->s2[idx] + 
		           gy_vuu->s2[idx] * gy_up->s2[idx] + gz_vuu->s2[idx] * gz_up->s2[idx] );
	    vxc_up[idx] -=  vxc2_upup[idx] * d2rho_up->s2[idx] ;
	    
	    vxc_up[idx] -= ( gx_vud->s2[idx] * gx_dw->s2[idx] + 
		           gy_vud->s2[idx] * gy_dw->s2[idx] + gz_vud->s2[idx] * gz_dw->s2[idx] );
	    vxc_up[idx] -= vxc2_updw[idx] * d2rho_dw->s2[idx];
    }

#endif    

    /* Release our memory */  
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
    

#endif
}                               /* end xcgga */
