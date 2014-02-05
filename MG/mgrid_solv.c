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
 *   void mgrid_solv(rmg_double_t *v_mat, rmg_double_t *f_mat, rmg_double_t *work, 
 *                   int dimx, int dimy, int dimz,
 *                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz,
 *                   int level, int *nb_ids, int max_levels, int *pre_cyc,
 *                   int *post_cyc, int mu_cyc, rmg_double_t step)
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
 *   k: When zero the call to solv_pois solves poissons equation. When non-zero helmholtz type
 * OUTPUT
 *   v_mat[(xdim+2)*(dimy+2)*(dimz+2)]: array to contain the solution
 * PARENTS
 *   get_vh.c mg_eig_state.c
 * CHILDREN
 *   trade_images.c solv_pois.c eval_residual.c mg_restrict.c mg_prolong.c solv_pois.c
 * SOURCE
 */


#include <stdio.h>
#include "const.h"
#include "grid.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "mg.h"
#include "main.h"


/*
 */



void mgrid_solv (rmg_double_t * v_mat, rmg_double_t * f_mat, rmg_double_t * work,
                 int dimx, int dimy, int dimz,
                 rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, rmg_double_t step, rmg_double_t k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim)
{
    int i;
    int cycl;
    int incx, incy, incz, size, idx;
    rmg_double_t scale;
    int ione = 1;
    int dx2, dy2, dz2, siz2;
    int ixoff, iyoff, izoff;
    rmg_double_t *resid, *newf, *newv, *newwork;


/* precalc some boundaries */
    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    resid = work + 2 * size;


    scale = 2.0 / (gridhx * gridhx * get_xside() * get_xside());
    scale = scale + (2.0 / (gridhy * gridhy * get_yside() * get_yside()));
    scale = scale + (2.0 / (gridhz * gridhz * get_zside() * get_zside()));
    scale = step / scale;
    trade_images (f_mat, dimx, dimy, dimz, nb_ids, FULL_FD);

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
        solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, k);


        /* trade boundary info */
        if ((level == max_levels) && (cycl == pre_cyc[level]-1)) {
            trade_images (v_mat, dimx, dimy, dimz, nb_ids, FULL_FD);
        }
        else {
            trade_images (v_mat, dimx, dimy, dimz, nb_ids, FULL_FD);
        }

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
    trade_images (resid, dimx, dimy, dimz, nb_ids, FULL_FD);


/* size for next smaller grid */
    dx2 = MG_SIZE (dimx, level, gxsize, gxoffset, pxdim, &ixoff, ct.boundaryflag);
    dy2 = MG_SIZE (dimy, level, gysize, gyoffset, pydim, &iyoff, ct.boundaryflag);
    dz2 = MG_SIZE (dimz, level, gzsize, gzoffset, pzdim, &izoff, ct.boundaryflag);

    siz2 = (dx2 + 2) * (dy2 + 2) * (dz2 + 2);

/* set storage pointers in the current workspace */
    newv = &work[0];
    newf = &work[siz2];
    newwork = &work[2 * siz2];


    for (i = 0; i < mu_cyc; i++)
    {

        mg_restrict (resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

        /* call mgrid solver on new level */
        mgrid_solv (newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                    gridhy * 2.0, gridhz * 2.0, level + 1, nb_ids,
                    max_levels, pre_cyc, post_cyc, mu_cyc, step, k,
                    gxsize, gysize, gzsize,
                    gxoffset, gyoffset, gzoffset,
                    pxdim, pydim, pzdim);


        mg_prolong (resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        scale = ONE;

        QMD_daxpy (size, scale, resid, ione, v_mat, ione);

        /* re-solve on this grid level */

        trade_images (v_mat, dimx, dimy, dimz, nb_ids, FULL_FD);

        for (cycl = 0; cycl < post_cyc[level]; cycl++)
        {

            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, k);

            /* trade boundary info */
            trade_images (v_mat, dimx, dimy, dimz, nb_ids, FULL_FD);

        }                       /* end for */

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
            trade_images (resid, dimx, dimy, dimz, nb_ids, FULL_FD);

        }                       /* end if */


    }                           /* for mu_cyc */

}

/******/
