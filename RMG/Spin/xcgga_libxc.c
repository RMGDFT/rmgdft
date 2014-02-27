/************************** SVN Revision Information **************************
 **    $Id$    **
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


#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>
#include "../../lib/libxc/lib64/include/xc.h"   
/* include Libxc's header file */

void xcgga_libxc (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int mode)
{

    int idx, sizr;
    rmg_double_t *gx, *gy, *gz, *vgx, *vgy, *vgz, *d2rho;
    int func_id_x, func_id_c;
    xc_func_type func_x, func_c;
    rmg_double_t *ec, *vc, *vsigma, *vsigma_c, *sigma;
    rmg_double_t hxxgrid, hyygrid, hzzgrid;

    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    FP0_BASIS = get_FP0_BASIS();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();

    sizr = FP0_BASIS;

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
        error_handler("Functional %d not found\n", func_id_x);
    
    if(xc_func_init(&func_c, func_id_c, XC_UNPOLARIZED) != 0)
        error_handler("Functional %d not found\n", func_id_c);

    my_calloc (ec, 5 * FP0_BASIS, rmg_double_t);
    vc = ec + FP0_BASIS; 
    sigma = ec + 2 * FP0_BASIS; 
    vsigma = ec + 3 * FP0_BASIS; 
    vsigma_c = ec + 4 * FP0_BASIS; 
    
    /* Grab some memory */
    my_malloc (gx, sizr, rmg_double_t);
    my_malloc (gy, sizr, rmg_double_t);
    my_malloc (gz, sizr, rmg_double_t);
    my_malloc (vgx, sizr, rmg_double_t);
    my_malloc (vgy, sizr, rmg_double_t);
    my_malloc (vgz, sizr, rmg_double_t);
//    my_malloc (agg, sizr, rmg_double_t);
    my_malloc (d2rho, sizr, rmg_double_t);




    /* Generate the gradient of the density */
    app_grad (rho, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);


    /* Get the Laplacian of the density */
    app6_del2 (rho, d2rho, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid );


    /* Absolute value of grad(rho) \dot grad(rho)*/
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        sigma[idx] =(gx[idx] * gx[idx] +
                           gy[idx] * gy[idx] + gz[idx] * gz[idx]);

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
    app_grad (vsigma, vgx, vgy, vgz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);


     /* Vxc = vrho -2 \div \dot ( vsigma * \grad(rho) ) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	     vxc[idx] -= 2.0 * ( vgx[idx] * gx[idx] + 
	     		   vgy[idx] * gy[idx] + vgz[idx] * gz[idx] ) ;
	     vxc[idx] -= 2.0 * vsigma[idx] * d2rho[idx];
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
