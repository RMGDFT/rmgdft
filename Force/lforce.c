/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



void lforce (double * rho, double * vh, double *force)
{

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    double r, Zv, rc, rc2, rcnorm, t1;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;

    double bx[3], norm1, fx, fy, fz;


    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;

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

    double *gx, *gy, *gz;
    gx = (double *)malloc(3*FP0_BASIS * sizeof(double));
    gy = gx + FP0_BASIS;
    gz = gx + 2*FP0_BASIS;

    int ithree = 3;
    double alpha = -get_vel_f(), zero = 0.0, mone = -1.0, *force_tmp;
    force_tmp = (double *)malloc(pct.num_loc_ions * 3 * sizeof(double));


    app_grad (vh, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);

    dgemm("T", "N", &ithree, &pct.num_loc_ions, &FP0_BASIS, &alpha, gx, &FP0_BASIS, 
            pct.localrhoc, &FP0_BASIS, &zero, force_tmp, &ithree); 


    app_grad (rho, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);

    dgemm("T", "N", &ithree, &pct.num_loc_ions, &FP0_BASIS, &alpha, gx, &FP0_BASIS, 
           pct.localpp, &FP0_BASIS, &mone, force_tmp, &ithree); 


    int ion1;
    for(ion1 = 0; ion1 <pct.num_loc_ions; ion1++)
    {
        ion = pct.loc_ions_list[ion1];
        
        force[ion *3 + 0] = force_tmp[ion1 *3 + 0];
        force[ion *3 + 1] = force_tmp[ion1 *3 + 1];
        force[ion *3 + 2] = force_tmp[ion1 *3 + 2];

    }
    

    free(force_tmp);
    free(gx);


}                               /* end lforce */

/******/
