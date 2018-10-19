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
 *   void mgrid_solv(double *v_mat, double *f_mat, double *Work, 
 *                   int dimx, int dimy, int dimz,
 *                   double gridhx, double gridhy, double gridhz,
 *                   int level, int *nb_ids, int max_levels, int *pre_cyc,
 *                   int *post_cyc, int mu_cyc, double step)
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


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "transition.h"
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void MgridSolvLocal(double * v_mat, double * f_mat, double * work,
                int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                int *post_cyc, double step, 
                int mu_cyc, int istate, int *iion, double Zfac)
{
    int i;
    int ione = 1;
    int dx2, dy2, dz2, siz2;
    double *resid, *newf, *newv, *newwork;
    int ib = 1;
//    if(level == 0) ib = 4;
//    if(level == 1) ib = 1;

    FiniteDiff *FD = new FiniteDiff(&Rmg_L);
//printf("LLLLLLLLL %d  %d\n",level,max_levels);
    int ixoff = 0;
    int iyoff = 0;
    int izoff = 0;
/* precalc some boundaries */
    int size = (dimx + 2) * (dimy + 2) * (dimz + 2);

    resid = work + 2 * size;

    double scale = 2.0 / (gridhx * gridhx * get_xside() * get_xside());
    scale = scale + (2.0 / (gridhy * gridhy * get_yside() * get_yside()));
    scale = scale + (2.0 / (gridhz * gridhz * get_zside() * get_zside()));
    scale = step / (scale + Zfac);

    for (int idx = 0; idx < size; idx++)
    {
        v_mat[idx] = 0.5 * scale * f_mat[idx];
    }
    ZeroBoundary(v_mat, dimx,dimy,dimz, ib);


/*
 * solve on this grid level 
 */

    for (int cycl = 0; cycl < pre_cyc[level]; cycl++)
    {
        /* solve once */
        SolvPoisLocal(FD, v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, 0);
    }



/*
 * on coarsest grid, we are finished
 */

    if (level >= max_levels)
    {
        delete FD;
        return;
    }


/* evaluate residual */
    EvalResidualLocal(FD, v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);

/* size for next smaller grid */
    dx2 = dimx / 2 + 1;
    dy2 = dimy / 2 + 1;
    dz2 = dimz / 2 + 1;
    siz2 = (dx2 + 2) * (dy2 + 2) * (dz2 + 2);
    if(!(dx2 % 2)) max_levels = 0;
    if(!(dy2 % 2)) max_levels = 0;
    if(!(dz2 % 2)) max_levels = 0;

/* set storage pointers in the current workspace */
    newv = &work[0];
    newf = &work[siz2];
    newwork = &work[2 * siz2];

    for (i = 0; i < mu_cyc; i++)
    {

        //ZeroBoundary(resid, dimx,dimy,dimz, ib);
        mg_restrict(resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        ////ZeroBoundary(newf, dx2,dy2,dz2, 1);

        /* call mgrid solver on new level */
        MgridSolvLocal(newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                   gridhy * 2.0, gridhz * 2.0, level + 1, nb_ids,
                   max_levels, pre_cyc, post_cyc, step,
                   1, istate, iion, 2.0*Zfac);

        ////ZeroBoundary(newv, dx2,dy2,dz2, 1);
        mg_prolong(resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        ZeroBoundary(resid, dimx,dimy,dimz, ib);


        scale = ONE;

        QMD_daxpy(size, scale, resid, ione, v_mat, ione);
        //ZeroBoundary(v_mat, dimx,dimy,dimz, ib);

        /* re-solve on this grid level */

        for (int cycl = 0; cycl < post_cyc[level]; cycl++)
        {
            /* solve once */
            SolvPoisLocal(FD, v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, 0);
            ZeroBoundary(v_mat, dimx,dimy,dimz, ib);
        }                       /* end for */

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            EvalResidualLocal(FD, v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
            ZeroBoundary(resid, dimx,dimy,dimz, ib);

        }                       /* end if */


    }                           /* for mu_cyc */

    delete FD;
}

void EvalResidualLocal (FiniteDiff *FD, double *mat, double *f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, double *res)
{
    int size = (dimx + 2) * (dimy + 2) * (dimz + 2);

    for (int idx = 0; idx < size; idx++) res[idx] = 0.0;

    FD->app2_del2 (mat, res, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    for (int idx = 0; idx < size; idx++) res[idx] = f_mat[idx] - res[idx];
}

void SolvPoisLocal (FiniteDiff *FD, double *vmat, double *fmat, double *work,
                int dimx, int dimy, int dimz, double gridhx,
                double gridhy, double gridhz, double step, double Zfac, double k)
{
    int size = (dimx + 2) * (dimy + 2) * (dimz + 2);

    for (int idx = 0; idx < size; idx++) work[idx] = 0.0;

    double diag = -FD->app2_del2(vmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    double scale = 1.0 / (diag + Zfac);
    scale = step * scale;

    // Non-zero k effectively means we are solving the Helmholtz rather than Poissons equation
    if(k != 0.0) {

        for (int idx = 0; idx < size; idx++)
        {
            vmat[idx] += scale * (work[idx] - k*vmat[idx] - fmat[idx]);
        }

     }
     else {

        for (int idx = 0; idx < size; idx++) vmat[idx] += scale * (work[idx] - fmat[idx]);

     }

}

