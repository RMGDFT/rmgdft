#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>



void get_vxc (double * rho, double * rho_oppo, double * rhocore, double * vxc)
{

    int idx;
    int FP0_BASIS;
    double *exc, *nrho, *nrho_oppo;

    FP0_BASIS = get_FP0_BASIS();
 
    my_malloc (exc, FP0_BASIS, double);
    my_malloc (nrho, FP0_BASIS, double);
    my_malloc (nrho_oppo, FP0_BASIS, double);


    /* Add in the nonlinear core charges correction from pseudopotential file */
    if (ct.spin_flag)
    {
        /*In spin polarized calculation,  add half of the nonlinear core charge for both 
	 * processor's own spin density and opposite spin density */
    	for (idx = 0; idx < FP0_BASIS; idx++)
	{
        	nrho[idx] = rhocore[idx] * 0.5 + rho[idx];         
		nrho_oppo[idx] = rhocore[idx] * 0.5 + rho_oppo[idx];
	}
    }
    else
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
        	nrho[idx] = rhocore[idx] + rho[idx];
    }

    get_vxc_exc(nrho, nrho_oppo, vxc, exc, ct.xctype);


    my_free (nrho_oppo);
    my_free (nrho);
    my_free (exc);


}                               /* end get_vxc */



/******/
