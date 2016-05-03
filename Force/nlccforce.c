/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/nlccforce.c *****
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
 *   void nlccforce(double *rho, double *vxc)
 *   Evaluates the ionic force component due to the non-linear core
 *   correction terms. 
 * INPUTS
 *   rho: total charge density
 *   vxc: exchange correlation potential 
 *   core charge density is obtained from ct.sp->rhocorelig
 * OUTPUT
 *   forces are added to structure ct.ion
 * PARENTS
 *   force.c
 * CHILDREN
 *   get_index.c to_crystal.c to_cartesian.c linint.c
 * SOURCE
 */




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



void nlccforce (double * rho, double * vxc, double *force_nlcc)
{

    double sumx, sumy, sumz;

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    double r;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;

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
    double rho_at_r;
    gx = (double *)malloc(3*FP0_BASIS * sizeof(double));
    gy = gx + FP0_BASIS;
    gz = gx + 2*FP0_BASIS;


    app_grad (vxc, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);

    int ithree = 3;
    double alpha = -get_vel_f(), zero = 0.0, *force_tmp;
    
    force_tmp = (double *)malloc(pct.num_loc_ions * 3 * sizeof(double));

    dgemm("T", "N", &ithree, &pct.num_loc_ions, &FP0_BASIS, &alpha, gx, &FP0_BASIS, 
            pct.localrhonlcc, &FP0_BASIS, &zero, force_tmp, &ithree); 
    int ion1;
    for(ion1 = 0; ion1 <pct.num_loc_ions; ion1++)
    {
        ion = pct.loc_ions_list[ion1];
        force_nlcc[ion *3 + 0] = force_tmp[ion1 *3 + 0];
        force_nlcc[ion *3 + 1] = force_tmp[ion1 *3 + 1];
        force_nlcc[ion *3 + 2] = force_tmp[ion1 *3 + 2];

    }
    

    free(force_tmp);
    free(gx);

}                               /* end nlccforce */

/******/
