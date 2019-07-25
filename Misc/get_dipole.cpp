/************************** SVN Revision Information **************************
  **    $Id: get_dipole.c 1643 2011-12-02 06:01:27Z miro $    **
******************************************************************************/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common_prototypes.h"
#include "main.h"

/*Conversion factor to convert dipole in a.u. to debye*/
#define DEBYE_CONVERSION 2.54174618772314


void get_dipole (double * rho)
{

    double px = 0.0, py = 0.0, pz = 0.0;
    double hxgrid, hygrid, hzgrid;
    double xc, yc, zc, ax[3], bx[3], x, y, z, temp, icharge, vel;
    int j, ix, iy, iz, ion, FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, dimx, dimy, dimz ;
    ION *iptr;


    vel = get_vel_f();

    hxgrid = get_hxxgrid();
    hygrid = get_hyygrid();
    hzgrid = get_hzzgrid();

    dimx = get_FPX0_GRID();
    dimy = get_FPY0_GRID();
    dimz = get_FPZ0_GRID();

    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();


    /*This gets dipole moment for electrons */
    xc = FPX_OFFSET * hxgrid;

    j = 0;
    for (ix = 0; ix < dimx; ix++)
    {

        yc = FPY_OFFSET * hygrid;
        for (iy = 0; iy < dimy; iy++)
        {

            zc = FPZ_OFFSET * hzgrid;
            for (iz = 0; iz < dimz; iz++)
            {

                /*vector between a grid point and a middle of the cell in crystal coordinates */
                ax[0] = xc - 0.5;
                ax[1] = yc - 0.5;
                ax[2] = zc - 0.5;

                /*Transform ax into vector in cartesian coordinates */
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                /*temp = rho[ix*incx + iy*incy +iz] * c->vel; */
                temp = rho[j];

                px += x * temp;
                py += y * temp;
                pz += z * temp;

                zc += hzgrid;
                j++;

            }                   /* end for */

            yc += hygrid;

        }                       /* end for */

        xc += hxgrid;

    }                           /* end for */





    /* Sum these up over all processors, multiply by volume elemnt so that sum is 
     * integration and invert sign, since electron charge should be negative*/
    px = -1.0 * vel * real_sum_all (px, pct.img_comm);
    py = -1.0 * vel * real_sum_all (py, pct.img_comm);
    pz = -1.0 * vel * real_sum_all (pz, pct.img_comm);


    /*Now we have dipole moment for electrons, need to add ions now */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Ionic pointer and ionic charge */
        iptr = &ct.ions[ion];
        icharge = ct.sp[iptr->species].zvalence;

        /*Difference vector between center of the cell and ionic position */
	ax[0] = iptr->xtal[0] - 0.5;
	ax[1] = iptr->xtal[1] - 0.5;
	ax[2] = iptr->xtal[2] - 0.5;



        /*Find vector ax in cartesian coordinates */
        to_cartesian (ax, bx);

        /*Ionic contribution to dipole moment */
        px += icharge * bx[0];
        py += icharge * bx[1];
        pz += icharge * bx[2];

    }                           /*end for (ion=0; ion<c->num_ions; ion++) */

    /*Now we need to convert to debye units */
    if (pct.imgpe==0)
    {
	printf("\n\n Dipole moment [Debye]: Absolute value: %.3f, vector: (%.3f,%.3f, %.3f)", 
		DEBYE_CONVERSION * sqrt (px * px + py * py + pz * pz), 
		DEBYE_CONVERSION *px, 
		DEBYE_CONVERSION *py, 
		DEBYE_CONVERSION *pz);
    }




}                               /* end getpoi_bc.c */

/******/
