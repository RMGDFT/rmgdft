#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>
//#include "../../lib/libxc/lib64/include/xc.h"   
#include "xc.h"   
/* include Libxc's header file */

void mgga_libxc (double * rho, double * tau, double * vxc, double * exc, int mode)
{

    int idx, sizr;
    double *gx, *gy, *gz, *vgx, *vgy, *vgz, *d2rho;
    int func_id_x, func_id_c;
    xc_func_type func_x, func_c;
    double *ec, *vc, *vsigma, *sigma;
    double hxxgrid, hyygrid, hzzgrid;
    
    double *vtau, *vlapl_rho;

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
    if (mode == MGGA_TB09)
    {
	    func_id_x = 208;   /* XC_mGGA_X_TB09 = 208 */
	    func_id_c = 9;   /* XC_LDA_C_PZ = 9 */
    }

    /* initialization */
    if(xc_func_init(&func_x, func_id_x, XC_UNPOLARIZED) != 0)
        error_handler("Functional %d not found\n", func_id_x);
    
    if(xc_func_init(&func_c, func_id_c, XC_UNPOLARIZED) != 0)
        error_handler("Functional %d not found\n", func_id_c);

    my_calloc (ec, 4 * FP0_BASIS, double);
    vc = ec + FP0_BASIS; 
    sigma = ec + 2 * FP0_BASIS; 
    vsigma = ec + 3 * FP0_BASIS; 


    my_calloc (vtau, 2 * FP0_BASIS, double);
    vlapl_rho = vtau + FP0_BASIS; 
    
    /* Grab some memory */
    my_malloc (gx, sizr, double);
    my_malloc (gy, sizr, double);
    my_malloc (gz, sizr, double);
    my_malloc (vgx, sizr, double);
    my_malloc (vgy, sizr, double);
    my_malloc (vgz, sizr, double);
//    my_malloc (agg, sizr, double);
    my_malloc (d2rho, sizr, double);




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


    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	if (tau[idx] <0)
		printf("idx: %d, rho: %16.9f, tau: %16.9e\n", idx, rho[idx], tau[idx]);
    } 





    /* get the exchange part for energy and potential first*/    
    if (func_x.info->family == XC_FAMILY_MGGA)
        xc_mgga_vxc (&func_x, FP0_BASIS, rho, sigma, d2rho, tau, vxc, vsigma, vlapl_rho, vtau); 


    /* get the correlation part for energy and potential */    
    if (func_c.info->family == XC_FAMILY_LDA)
        xc_lda_exc_vxc (&func_c, FP0_BASIS, rho, ec, vc); 
    
    xc_func_end (&func_x);
    xc_func_end (&func_c); 

    
    /* add exchange correlation together */ 
    for (idx = 0; idx < FP0_BASIS; idx++) 
    {
	    vxc[idx] += vc[idx];
	    exc[idx] += ec[idx];
    }




    my_free (ec);
    
    my_free (vtau);

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
