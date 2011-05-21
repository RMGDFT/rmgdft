/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/xcgga_spin_libxc.c *****
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

#define SMALL 1.e-8


void xcgga_spin_libxc(REAL * rho_up, REAL * rho_dw, REAL * vxc_up, REAL * exc, int mode) 
{

    int idx;
    FP0_GRID *d2rho_up, *d2rho_dw;
    FP0_GRID *gx_up, *gy_up, *gz_up;
    FP0_GRID *gx_dw, *gy_dw, *gz_dw; 
    FP0_GRID *gx_vsigmauu, *gy_vsigmauu, *gz_vsigmauu, *gx_vsigmaud, *gy_vsigmaud, *gz_vsigmaud; 
    REAL rhospin[FP0_BASIS][2], sigma[FP0_BASIS][3], vsigma[FP0_BASIS][3], vspin[FP0_BASIS][2];
    int func_id_x, func_id_c;
    xc_func_type func_x, func_c;
    REAL *ec, *vsigma_upup, *vsigma_updw;

    
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
    if(xc_func_init(&func_x, func_id_x, XC_POLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_x);
    	exit(1);
    }
    
    if(xc_func_init(&func_c, func_id_c, XC_POLARIZED) != 0)
    {
	if (pct.imgpe == 0)    
    		printf( "Functional %d not found\n", func_id_c);
    	exit(1);
    } 


    my_calloc (ec, 3 * FP0_BASIS, REAL);
    vsigma_upup = ec + FP0_BASIS; 
    vsigma_updw = ec + 2 * FP0_BASIS; 


    /* Grab some memory */ 
    /* to hold gradient of spin up charge density */ 
    my_malloc (gx_up, 1, FP0_GRID);
    my_malloc (gy_up, 1, FP0_GRID);
    my_malloc (gz_up, 1, FP0_GRID);

    /* to hold gradient of spin down charge density */
    my_malloc (gx_dw, 1, FP0_GRID);
    my_malloc (gy_dw, 1, FP0_GRID);
    my_malloc (gz_dw, 1, FP0_GRID);



    
    /* to hold laplaciant of the spin up and down charge density */
    my_malloc (d2rho_up, 1, FP0_GRID);
    my_malloc (d2rho_dw, 1, FP0_GRID); 

    /* to hold the gradient of potentials */
    my_malloc (gx_vsigmauu, 1, FP0_GRID);
    my_malloc (gy_vsigmauu, 1, FP0_GRID);
    my_malloc (gz_vsigmauu, 1, FP0_GRID);

    my_malloc (gx_vsigmaud, 1, FP0_GRID);
    my_malloc (gy_vsigmaud, 1, FP0_GRID);
    my_malloc (gz_vsigmaud, 1, FP0_GRID);


    /* Generate the gradient of the density */
    app_gradf (rho_up, gx_up, gy_up, gz_up);    /* spin up density */
    app_gradf (rho_dw, gx_dw, gy_dw, gz_dw);    /* spin down density */
    

    /* Get the Laplacian of the density */
    app6_del2f (rho_up, d2rho_up);
    app6_del2f (rho_dw, d2rho_dw);

    
    /* pack rho and rho_oppo into the 2D array rhospin */ 
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    rhospin[idx][0] = rho_up[idx];
	    rhospin[idx][1] = rho_dw[idx];
    }


    /* pack the 2D array sigma[][2] */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        
	/* grad(rho_up) \dot grad(rho_up) */    
        sigma[idx][0] =  (gx_up->s2[idx] * gx_up->s2[idx] +
                             gy_up->s2[idx] * gy_up->s2[idx] + gz_up->s2[idx] * gz_up->s2[idx]);
        
	/* grad(rho_up) \dot grad(rho_dw)*/
        sigma[idx][1] =  (gx_up->s2[idx] * gx_dw->s2[idx] +
                             gy_up->s2[idx] * gy_dw->s2[idx] + gz_up->s2[idx] * gz_dw->s2[idx]);
        
	/*  (grad rho_dw) \dot (grad rho_dw)  */
        sigma[idx][2] = gx_dw->s2[idx] * gx_dw->s2[idx] +
                	             gy_dw->s2[idx] * gy_dw->s2[idx] + gz_dw->s2[idx] * gz_dw->s2[idx];

    }                           
   
    /* get the exchange part for energy and potential first*/    
    if (func_x.info->family == XC_FAMILY_GGA)
        xc_gga_exc_vxc (&func_x, FP0_BASIS, &rhospin[0][0], &sigma[0][0], exc, &vspin[0][0], &vsigma[0][0]); 




    /* unpack the 2D array of exchange potential to vxc_up, vsigma_upup, vsigma_updw */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	    vxc_up[idx] = vspin[idx][0];
	    vsigma_upup[idx] = vsigma[idx][0];
	    vsigma_updw[idx] = vsigma[idx][1];
    }

    /* get the correlation part for energy and potential */    
    if (func_c.info->family == XC_FAMILY_GGA)
        xc_gga_exc_vxc (&func_c, FP0_BASIS, &rhospin[0][0], &sigma[0][0], ec, &vspin[0][0], &vsigma[0][0]); 
 

    
    /* add up the correlation part together */
    for (idx = 0; idx < FP0_BASIS; idx++) 
    {
	    vxc_up[idx] += vspin[idx][0];
	    exc[idx] += ec[idx];
	    vsigma_upup[idx] += vsigma[idx][0];
	    vsigma_updw[idx] += vsigma[idx][1];
    }

    xc_func_end (&func_x);
    xc_func_end (&func_c); 

    /* Get gradient of vsigma_upup and vsigma_updw */
    app_gradf (vsigma_upup, gx_vsigmauu, gy_vsigmauu, gz_vsigmauu); 
    app_gradf (vsigma_updw, gx_vsigmaud, gy_vsigmaud, gz_vsigmaud); 




     /* Vxc_up = vrho_up - \div \dot ( 2* vsigma_upup * \grad(rho_up) + vsigma_updw * \grad(rho_dw) ) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	     vxc_up[idx] -= 2.0 * ( gx_vsigmauu->s2[idx] * gx_up->s2[idx] + 
	     		   gy_vsigmauu->s2[idx] * gy_up->s2[idx] + gz_vsigmauu->s2[idx] * gz_up->s2[idx] ) ;
	     vxc_up[idx] -= 2.0 * vsigma_upup[idx] * d2rho_up->s2[idx];
	     vxc_up[idx] -=  ( gx_vsigmaud->s2[idx] * gx_dw->s2[idx] + 
	     		   gy_vsigmaud->s2[idx] * gy_dw->s2[idx] + gz_vsigmaud->s2[idx] * gz_dw->s2[idx] ) ;
	     vxc_up[idx] -=  vsigma_updw[idx] * d2rho_dw->s2[idx]; 
    }




    my_free (ec);

    /* Release our memory */  
    my_free (gx_up);
    my_free (gy_up);
    my_free (gz_up);

    my_free (gx_dw);
    my_free (gy_dw);
    my_free (gz_dw);
    
    my_free(d2rho_up);
    my_free(d2rho_dw);

    my_free (gz_vsigmauu);
    my_free (gy_vsigmauu);
    my_free (gx_vsigmauu);

    my_free (gz_vsigmaud);
    my_free (gy_vsigmaud);
    my_free (gx_vsigmaud);


}
