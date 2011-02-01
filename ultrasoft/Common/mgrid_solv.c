/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/mgrid_solv.c *****
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
 *   void mgrid_solv(REAL *v_mat, REAL *f_mat, REAL *work, 
 *                   int dimx, int dimy, int dimz,
 *                   REAL gridhx, REAL gridhy, REAL gridhz,
 *                   int level, int *nb_ids, int max_levels, int *pre_cyc,
 *                   int *post_cyc, int mu_cyc, REAL step)
 *	solves the Poisson equation: del^2 v = f for v. 
 *      This routine is called recursively for each level of multi-grid.
 * INPUTS
 *   f_mat[(xdim+2)*(dimy+2)*(dimz+2)]: the charge density (or residual) data
 *   work: workspace, dimensioned at least as large as the above matricies
 *   dimx,dimy,dimz:  size of the data arrays
 *   gridhx:  grid spacing for the x plane
 *   gridhy:  grid spacing for the y plane
 *   gridhz:  grid spacing for the z plane
 *   level:   indicates current level of multigrid
 *   nb_ids:  PE neighbor list
 *   max_levels:  maximum multigrid level
 *   pre_cyc: 
 *   post_cyc: number of post-smoothing
 *   mu_cyc:   number of iteration on the same level
 *   step:  time step
 * OUTPUT
 *   v_mat[(xdim+2)*(dimy+2)*(dimz+2)]: array to contain the solution
 * PARENTS
 *   get_vh.c mg_eig_state.c
 * CHILDREN
 *   trade_images.c solv_pois.c eval_residual.c mg_restrict.c mg_prolong.c solv_pois.c
 * SOURCE
 */


#include <stdio.h>

#include "main.h"

/*
 */



void mgrid_solv (REAL * v_mat, REAL * f_mat, REAL * work,
                 int dimx, int dimy, int dimz,
                 REAL gridhx, REAL gridhy, REAL gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, REAL step)
{
    int i;
    int cycl;
    int incx, incy, incz, size, idx;
    REAL scale;
    int ione = 1;
    int dx2, dy2, dz2, siz2;
    REAL *resid, *newf, *newv, *newwork;


/* precalc some boundaries */
    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    resid = work + 2 * size;


    scale = 2.0 / (gridhx * gridhx * ct.xside * ct.xside);
    scale = scale + (2.0 / (gridhy * gridhy * ct.yside * ct.yside));
    scale = scale + (2.0 / (gridhz * gridhz * ct.zside * ct.zside));
    scale = step / scale;
    trade_images (f_mat, dimx, dimy, dimz, nb_ids);


    for (idx = 0; idx < size; idx++)
    {

        v_mat[idx] = -scale * f_mat[idx];

    }                           /* end for */


/*
 * solve on this grid level 
 */


    for (cycl = 0; cycl < pre_cyc[level]; cycl++)
    {

        /* solve once */
        solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step);


        /* trade boundary info */
        trade_images (v_mat, dimx, dimy, dimz, nb_ids);

    }


/*
 * on coarsest grid, we are finished
 */

    if (level == max_levels)
    {

        return;

    }                           /* end if */


/* evaluate residual */
    eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
    trade_images (resid, dimx, dimy, dimz, nb_ids);


/* size for next smaller grid */
    dx2 = dimx / 2;
    dy2 = dimy / 2;
    dz2 = dimz / 2;
    siz2 = (dx2 + 2) * (dy2 + 2) * (dz2 + 2);

/* set storage pointers in the current workspace */
    newv = &work[0];
    newf = &work[siz2];
    newwork = &work[2 * siz2];


    for (i = 0; i < mu_cyc; i++)
    {

        mg_restrict (resid, newf, dimx, dimy, dimz);

        /* call mgrid solver on new level */
        mgrid_solv (newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                    gridhy * 2.0, gridhz * 2.0, level + 1, nb_ids,
                    max_levels, pre_cyc, post_cyc, 1, step);


        mg_prolong (resid, newv, dimx, dimy, dimz);
        scale = ONE;

        QMD_saxpy (size, scale, resid, ione, v_mat, ione);

        /* re-solve on this grid level */

        trade_images (v_mat, dimx, dimy, dimz, nb_ids);

        for (cycl = 0; cycl < post_cyc[level]; cycl++)
        {

            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step);

            /* trade boundary info */
            trade_images (v_mat, dimx, dimy, dimz, nb_ids);

        }                       /* end for */

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
            trade_images (resid, dimx, dimy, dimz, nb_ids);

        }                       /* end if */


    }                           /* for mu_cyc */

}

/******/
