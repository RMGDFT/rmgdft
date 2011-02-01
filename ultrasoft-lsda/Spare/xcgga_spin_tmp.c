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


#define    crs    1.91915829267751281
//#define    SMALL  1.e-8
#define    SMALL  1.e-18

/* up and dn are convenient for naming the processor's own spin the the opposite spin, 
 doesn't have real physical meaning*/

void xcgga_spin_tmp(REAL * rho_up, REAL * rho_dw, REAL * vxc_up, REAL * exc, int mode)

{
#if 1
    int ix, iy, iz, idx;
    FP0_GRID *gx, *gy, *gz, *agx, *agy, *agz, *agg, *d2rho, *d2rho_up, *d2rho_dw;
    FP0_GRID *gx_up, *gy_up, *gz_up, *agg_up, *agx_up, *agy_up, *agz_up;
    FP0_GRID *gx_dw, *gy_dw, *gz_dw, *agg_dw, *agx_dw, *agy_dw, *agz_dw;
    REAL d, kf, us, uu, t, vv, ww;
    REAL gdgz, rho2, s_up, us_up, v_up, d_up, u_up, s_dw, us_dw, v_dw, d_dw, u_dw, kf_up, kf_dw;

    REAL pisq3, ex, vx, ec;
    REAL zet, rs, g, h, sk, gks2;
    REAL dvc_up, dvc_dw;
    REAL ecrs, eczet, alfc;
    REAL cen, xen, cpot;
    REAL cpot_up, cpot_dw;
   
    REAL vx_up, vx_dw, vc_up, vc_dw, ex_up, ex_dw;

    int ndim = 3, lgga, lpot;

    pisq3 = THREE * PI * PI;
    REAL rho_tot[FP0_BASIS], vxc_dw[FP0_BASIS];

    /* Grab some memory */
    my_malloc (gx, 1, FP0_GRID);
    my_malloc (gy, 1, FP0_GRID);
    my_malloc (gz, 1, FP0_GRID);
    my_malloc (agx, 1, FP0_GRID);
    my_malloc (agy, 1, FP0_GRID);
    my_malloc (agz, 1, FP0_GRID);
    my_malloc (agg, 1, FP0_GRID);
    my_malloc (d2rho, 1, FP0_GRID);

    my_malloc (agx_up, 1, FP0_GRID);
    my_malloc (agy_up, 1, FP0_GRID);
    my_malloc (agz_up, 1, FP0_GRID);
    my_malloc (agx_dw, 1, FP0_GRID);
    my_malloc (agy_dw, 1, FP0_GRID);
    my_malloc (agz_dw, 1, FP0_GRID);
    my_malloc (gx_up, 1, FP0_GRID);
    my_malloc (gy_up, 1, FP0_GRID);
    my_malloc (gz_up, 1, FP0_GRID);
    my_malloc (agg_up, 1, FP0_GRID);
    my_malloc (gx_dw, 1, FP0_GRID);
    my_malloc (gy_dw, 1, FP0_GRID);
    my_malloc (gz_dw, 1, FP0_GRID);
    my_malloc (agg_dw, 1, FP0_GRID);
    my_malloc (d2rho_up, 1, FP0_GRID);
    my_malloc (d2rho_dw, 1, FP0_GRID);


    /*Total charge density*/
    for (idx = 0; idx < FP0_BASIS; idx++)
	rho_tot[idx] = rho_up[idx] + rho_dw[idx]; 

    /* Generate the gradient of the density */
    app_gradf (rho_tot, gx, gy, gz);                /*total density*/
    app_gradf (rho_up, gx_up, gy_up, gz_up);    /*spin up density*/
    app_gradf (rho_dw, gx_dw, gy_dw, gz_dw);    /*spin down density*/
    

    /* Get the Laplacian of the density */
    app6_del2f (rho_tot, d2rho);
    app6_del2f (rho_up, d2rho_up);
    app6_del2f (rho_dw, d2rho_dw);

    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        agg->s2[idx] = sqrt (gx->s2[idx] * gx->s2[idx] +
                             gy->s2[idx] * gy->s2[idx] + gz->s2[idx] * gz->s2[idx]);
        
	/*Absolute value of grad(rho_up)*/    
        agg_up->s2[idx] = sqrt (gx_up->s2[idx] * gx_up->s2[idx] +
                             gy_up->s2[idx] * gy_up->s2[idx] + gz_up->s2[idx] * gz_up->s2[idx]);
        
	/*Absolute value of grad(rho_dw) */
        agg_dw->s2[idx] = sqrt (gx_dw->s2[idx] * gx_dw->s2[idx] +
                             gy_dw->s2[idx] * gy_dw->s2[idx] + gz_dw->s2[idx] * gz_dw->s2[idx]);
    }                           /* end for */

    /* Get its gradient */
    app_gradf (agg->s2, agx, agy, agz);
    app_gradf (agg_up->s2, agx_up, agy_up, agz_up);
    app_gradf (agg_dw->s2, agx_dw, agy_dw, agz_dw);

    /* Now get the potential */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        /* for exchange part */

        if (mode == GGA_PBE)
        {
            /* exchange potential from Perdew, Burke, Ernzerhof */
            lpot = 1;
            lgga = 1;


            /* do spin up part first */
            
	    rho2 = TWO * rho_up[idx];

	    if ( rho2 > SMALL)
	    {
                    kf_up = pow (pisq3 * rho2, 0.333333333333333);        
	   	    s_up = TWO * agg_up->s2[idx] / (TWO * kf_up * rho2);
               	    us_up = FOUR * ( gx_up->s2[idx] * agx_up->s2[idx] + gy_up->s2[idx] * agy_up->s2[idx] + gz_up->s2[idx] * agz_up->s2[idx]);
           	    u_up = us_up / (rho2 * rho2 * EIGHT * kf_up * kf_up * kf_up);
           	    v_up = TWO * d2rho_up->s2[idx] / (FOUR * rho2 * kf_up * kf_up);
           	    exchpbe (&rho2, &s_up, &u_up, &v_up, &lgga, &lpot, &ex_up, &vx_up);
	    }
	    else
	    {
		    ex_up = 0.0;
		    vx_up = 0.0;
	    }

            /* Now do the spin down part */   
            
	    rho2 = TWO * rho_dw[idx];

	    if ( rho2 > SMALL)
	    {
                    kf_dw = pow (pisq3 * rho2, 0.333333333333333);        
           	    s_dw = TWO * agg_dw->s2[idx] / (TWO * kf_dw * rho2);
           	    us_dw = FOUR * ( gx_dw->s2[idx] * agx_dw->s2[idx] + gy_dw->s2[idx] * agy_dw->s2[idx] + gz_dw->s2[idx] * agz_dw->s2[idx] );
           	    u_dw = us_dw / (rho2 * rho2 * EIGHT * kf_dw * kf_dw * kf_dw);
           	    v_dw = TWO * d2rho_dw->s2[idx] / (FOUR * rho2 * kf_dw * kf_dw);
           	    exchpbe (&rho2, &s_dw, &u_dw, &v_dw, &lgga, &lpot, &ex_dw, &vx_dw);
	    }
	    else 
	    {
		    ex_dw = 0.0;
		    vx_dw = 0.0;
	    }
	    
            /* get the exchage energy */   

	    if ( rho_tot[idx] > SMALL) 
               	    ex = (ex_up * rho_up[idx] + ex_dw * rho_dw[idx]) / rho_tot[idx];
	    else
		    ex = 0.0;

        }                       /* end if */


	/* for correlation part */

        if (mode == GGA_PBE)
        {     
	    d =rho_tot[idx];	
	    if (d <  SMALL )
	    {
		vxc_up[idx] = 0.0;
		exc[idx] = 0.0;
	    }
	    else
	    {
                kf = pow (pisq3 * d, 0.333333333333333); 

        	us = gx->s2[idx] * agx->s2[idx] + gy->s2[idx] * agy->s2[idx] + gz->s2[idx] * agz->s2[idx];
            
		/* relative polarization for spin polarized */
            	zet = (rho_up[idx]-rho_dw[idx]) / d;
	    
                if ( fabs(zet) > 1.0)
		{
			
            		vxc_up[idx] =  vx_up ;
            		vxc_dw[idx] =  vx_dw ;
            		exc[idx] =   ex;
			continue;
		}
		
              	rs = crs / kf;
            	sk = TWO * sqrt (kf / PI);
         
            	/* for spin polarized calculation*/
            	g = pow (1.0 + zet, 0.6666666666666666);
            	g += pow (1.0 - zet, 0.6666666666666666);
            	g = g / TWO;
            	gks2 = TWO * sk * g;
            	t = agg->s2[idx] / (d * gks2);
            	uu = us / (d * d * gks2 * gks2 * gks2);
            	vv = d2rho->s2[idx] / (d * gks2 * gks2);
	    
	    	/* calculate |grad up|^2 - |grad dw|^2 - zet * |grad rho|^2 */
	    	gdgz = - zet * agg->s2[idx] * agg->s2[idx];
	    	gdgz += agg_up->s2[idx] * agg_up->s2[idx] - agg_dw->s2[idx] * agg_dw->s2[idx];
	    
	    	/* ww for spin polarized calculation*/
            	ww = gdgz /(d * d * gks2 * gks2);
            	lpot = 1;
            	lgga = 1;
            	corpbe (&rs, &zet, &t, &uu, &vv, &ww, &lgga, &lpot, &ec, &vc_up, &vc_dw,
                    &h, &dvc_up, &dvc_dw);

           	/* corelation potential for spin up and down respectively*/
            	cpot_up = vc_up + dvc_up;
            	cpot_dw = vc_dw + dvc_dw;           
            	cen = ec + h;
         
	    
#if 0           
	   if ( (vxc_up[idx] > 0 && rho_up[idx] < 1.e-12) || (vxc_dw[idx] > 0 && rho_dw[idx] < 1.e-12))
	   {
		   pz_spin(rs, zet, &cen, &cpot_up, &cpot_dw);
	   }
#endif 
	    
            	/* get exchange-correlation potential for spin up and down*/
            	vxc_up[idx] =  ( cpot_up + vx_up );
            	vxc_dw[idx] =  ( cpot_dw + vx_dw );
            	exc[idx] =  ( cen + ex );
	    }
	    
        }

    }                           /* end for */

    /* Release our memory */
    my_free (d2rho);
    my_free (agg);
    my_free (agz);
    my_free (agy);
    my_free (agx);
    my_free (gz);
    my_free (gy);
    my_free (gx);

    my_free (agz_up);
    my_free (agy_up);
    my_free (agx_up);
    my_free (agz_dw);
    my_free (agy_dw);
    my_free (agx_dw);
    my_free (gx_up);
    my_free (gy_up);
    my_free (gz_up);
    my_free (agg_up);
    my_free (gx_dw);
    my_free (gy_dw);
    my_free (gz_dw);
    my_free (agg_dw);
    my_free(d2rho_up);
    my_free(d2rho_dw);
#endif
}                               /* end xcgga */
