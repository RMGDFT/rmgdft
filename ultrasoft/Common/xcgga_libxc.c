/************************** SVN Revision Information **************************
 **    $Id: xcgga.c 1288 2011-03-15 00:03:28Z miro $    **
******************************************************************************/

/****f* QMD-MGDFT/xcgga_libxc.c *****
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
#include "../../lib/libxc/include/xc.h"   
/* include Libxc's header file */

void xcgga_libxc (REAL * rho, REAL * vxc, REAL * exc, int mode)
{

    int ix, iy, iz, idx;
    FP0_GRID *gx, *gy, *gz, *vgx, *vgy, *vgz, *d2rho;
    int func_id_x, func_id_c;
    xc_func_type func_x, func_c;
    REAL *ec, *vc, *vsigma, *vsigma_c, *sigma;

    /* get the correct functional id for each type */
    if (mode == GGA_PBE)
    {
	    func_id_x = 101;   /* XC_GGA_X_PBE = 101 */
	    func_id_c = 130;   /* XC_GGA_C_PBE = 130 */
    }
    else if (mode == GGA_XB_CP)
    {
	    func_id_x = 106;   /* XC_GGA_X_B88 = 106, XC_GGA_X_B86 = 103, XC_GGA_X_B86_MGC = 105 */
	    func_id_c = 134;   /* XC_GGA_C_P86 = 132  if doesn't work set value equal 134 */
    }
    else if (mode == GGA_XP_CP)
    {
	    func_id_x = 109;   /* XC_GGA_X_PW91 = 109 */
	    func_id_c = 134;   /* XC_GGA_C_PW91 = 134 */
    }
    else if (mode == GGA_BLYP)
    {
	    func_id_x = 106;   /* XC_GGA_X_B88 = 106 */
	    func_id_c = 131;   /* XC_GGA_C_LYP = 131 */
    }

    /* initialization */
    if(xc_func_init(&func_x, func_id_x, XC_UNPOLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_x);
    	exit(1);
    }
    
    if(xc_func_init(&func_c, func_id_c, XC_UNPOLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_c);
    	exit(1);
    } 


    my_calloc (ec, 5 * FP0_BASIS, REAL);
    vc = ec + FP0_BASIS; 
    sigma = ec + 2 * FP0_BASIS; 
    vsigma = ec + 3 * FP0_BASIS; 
    vsigma_c = ec + 4 * FP0_BASIS; 
    

    /* Grab some memory */
    my_malloc (gx, 1, FP0_GRID);
    my_malloc (gy, 1, FP0_GRID);
    my_malloc (gz, 1, FP0_GRID);
    my_malloc (vgx, 1, FP0_GRID);
    my_malloc (vgy, 1, FP0_GRID);
    my_malloc (vgz, 1, FP0_GRID);
    /* my_malloc (agg, 1, FP0_GRID); */
    my_malloc (d2rho, 1, FP0_GRID); 



    /* Generate the gradient of the density */
    app_gradf (rho, gx, gy, gz);


    /* Get the Laplacian of the density */
    app6_del2f (rho, d2rho);


    /* Absolute value of grad(rho) \dot grad(rho)*/
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        sigma[idx] =(gx->s2[idx] * gx->s2[idx] +
                           gy->s2[idx] * gy->s2[idx] + gz->s2[idx] * gz->s2[idx]);

    } 



    /* get the exchange part for energy and potential first*/    
    if (func_x.info->family == XC_FAMILY_GGA)
        xc_gga_exc_vxc (&func_x, FP0_BASIS, rho, sigma, exc, vxc, vsigma); 


    /* get the correlation part for energy and potential */    
    if (func_c.info->family == XC_FAMILY_GGA)
        xc_gga_exc_vxc (&func_c, FP0_BASIS, rho, sigma, ec, vc, vsigma_c); 
    
    xc_func_end (&func_x);
    xc_func_end (&func_c); 

    
    /* add exchange correlation together */ 
    for (idx = 0; idx < FP0_BASIS; idx++) 
    {
	    vxc[idx] += vc[idx];
	    exc[idx] += ec[idx];
	    vsigma[idx] += vsigma_c[idx];
    }


    /* Get gradient of vsigma */
    app_gradf (vsigma, vgx, vgy, vgz);


     /* Vxc = vrho -2 \div \dot ( vsigma * \grad(rho) ) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	     vxc[idx] -= 2.0 * ( vgx->s2[idx] * gx->s2[idx] + 
	     		   vgy->s2[idx] * gy->s2[idx] + vgz->s2[idx] * gz->s2[idx] ) ;
	     vxc[idx] -= 2.0 * vsigma[idx] * d2rho->s2[idx];
    }


    my_free (ec);

    /* Release our memory */
    my_free (d2rho);
    /* my_free (agg); */
    my_free (vgz);
    my_free (vgy);
    my_free (vgx);
    my_free (gz);
    my_free (gy);
    my_free (gx);

}
