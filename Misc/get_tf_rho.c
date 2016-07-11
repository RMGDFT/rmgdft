/************************** SVN Revision Information **************************
 **    $Id: init_nuc.c 2975 2015-03-20 01:42:51Z ebriggs $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "rmg_error.h"
#include "common_prototypes.h"
#include "main.h"
#include "AtomicInterpolate.h"


void get_tf_rho (double * tf_rho)
{

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FP0_BASIS;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;
    int tf_num_ions;

    double r, t1;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    TF_ION *tf_iptr;
    int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
    double lradius = 9.0;
    double alpha, alpha0, q0, q, prefactor, prefactor0;


    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    FP0_BASIS = get_FP0_BASIS();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();
    FNX_GRID = get_FNX_GRID();
    FNY_GRID = get_FNY_GRID();
    FNZ_GRID = get_FNZ_GRID();


    ilow = FPX_OFFSET;
    jlow = FPY_OFFSET;
    klow = FPZ_OFFSET;
    ihi = ilow + FPX0_GRID;
    jhi = jlow + FPY0_GRID;
    khi = klow + FPZ0_GRID;

    tf_num_ions = ct.num_tfions;


#if 1


    /* Initialize the TF charge array */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	tf_rho[idx] = 0.0;
    }
#endif

    /* Loop over ions */
    for (ion = 0; ion < tf_num_ions; ion++)
    {
	tf_iptr = &ct.tf_ions[ion];
	to_crystal (tf_iptr->xtal, tf_iptr->crds);
	

#if 1
	/*Coefficients for gaussians representing charge density and ions of water molecules */
	q = tf_iptr->q;
	q0 = tf_iptr->q0;
	alpha = tf_iptr->alpha;
	alpha0 = tf_iptr->alpha0;

	//printf("\n TF ion %d: q:%f q0:%f alpha:%f alpha0:%f", ion, q, q0, alpha, alpha0);


	dimx =  lradius/(hxxgrid*xside);
	dimy =  lradius/(hyygrid*yside);
	dimz =  lradius/(hzzgrid*zside);

	dimx = dimx * 2 + 1;
	dimy = dimy * 2 + 1;
	dimz = dimz * 2 + 1;


	xstart =  tf_iptr->xtal[0] / hxxgrid - dimx/2;
	xend = xstart + dimx;
	ystart =  tf_iptr->xtal[1] / hyygrid - dimy/2;
	yend = ystart + dimy;
	zstart =  tf_iptr->xtal[2] / hzzgrid - dimz/2;
	zend = zstart + dimz;
	

	prefactor = q * pow( (alpha/3.14159265), 1.5);
	prefactor0 = q0 * pow( (alpha0/3.14159265), 1.5);


	for (ix = xstart; ix < xend; ix++)
	{
	    // fold the grid into the unit cell
	    ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
	    if(ixx >= ilow && ixx < ihi)
	    {

		for (iy = ystart; iy < yend; iy++)
		{
		    // fold the grid into the unit cell
		    iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
		    if(iyy >= jlow && iyy < jhi)
		    {
			for (iz = zstart; iz < zend; iz++)
			{
			    // fold the grid into the unit cell
			    izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
			    if(izz >= klow && izz < khi)
			    {

				idx = (ixx-ilow) * FPY0_GRID * FPZ0_GRID + (iyy-jlow) * FPZ0_GRID + izz-klow;

				/*if ((idx < 0) || (idx >= FP0_BASIS))
				{
				    printf("\n Warning: idx out pf bounds, should be 0-%d, but is %d", FP0_BASIS, idx);
				    rmg_error_handler ("idx out of bounds");
				}*/
				    
				x[0] = ix * hxxgrid -  tf_iptr->xtal[0];
				x[1] = iy * hyygrid -  tf_iptr->xtal[1];
				x[2] = iz * hzzgrid -  tf_iptr->xtal[2];
				r = metric (x);

				tf_rho[idx] += prefactor * exp (-1.0*alpha*r*r);
				tf_rho[idx] -= prefactor0 * exp (-1.0*alpha0*r*r);



			    }                           /* end for */

			}
		    }
		}
	    }

	}
#endif
    }
    
    
#if 1
    /* Set net TF charge to 0 */
    t1 = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	t1 += tf_rho[idx];
    }


    t1 = real_sum_all (t1, pct.grid_comm);  /* sum over pct.grid_comm  */
    
    
    if (pct.imgpe==0)
	printf("\nTotal TF charge is %.8e \n", t1*get_vel_f());


    t1 /= (FP0_BASIS*npes);
    
    for (idx = 0; idx < FP0_BASIS; idx++)
	tf_rho[idx] -= t1;


    /*Check again, disable if works correctly*/
    t1 = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
	t1 += tf_rho[idx];
    
    t1 = real_sum_all (t1, pct.grid_comm);  /* sum over pct.grid_comm  */
    
    if (pct.imgpe==0)
	printf("\nTotal TF charge after adjustment is %.8e\n", t1*get_vel_f());


#endif





}                               /* end init_nuc */

/******/
