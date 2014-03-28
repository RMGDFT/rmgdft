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
 *   void mgrid_solv(RmgType *v_mat, RmgType *f_mat, RmgType *work, 
 *                   int dimx, int dimy, int dimz,
 *                   RmgType gridhx, RmgType gridhy, RmgType gridhz,
 *                   int level, int *nb_ids, int max_levels, int *pre_cyc,
 *                   int *post_cyc, int mu_cyc, RmgType step)
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


#include "Mgrid.h"
#include "FiniteDiff.h"
#include "BlasWrappers.h"
#include "TradeImages.h"
#include "RmgTimer.h"

using namespace std;

template void Mgrid::mgrid_solv<float>(float*, float*, float*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<double>(double*, double*, double*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, int, int, int, int, int, int, int, int, int, int);


Mgrid::Mgrid(Lattice *lptr, TradeImages *tptr)
{
    L = lptr;
    T = tptr;
    level_flag = 0;
}

Mgrid::~Mgrid(void)
{
    if(level_flag)
        cout << "Warning: too many multigrid levels were requested " << level_flag << " times." << endl;
}

template <typename RmgType>
void Mgrid::mgrid_solv (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    RmgTimer RT0("Mgrid_solv");
    std::string timername = "Mgrid_solv: level " + std::to_string(level);
    RmgTimer RT1(timername.c_str());

    int i;
    int cycl;
    int size, idx;
    RmgType scale;
    int ione = 1;
    int dx2, dy2, dz2, siz2;
    int ixoff, iyoff, izoff;
    RmgType *resid, *newf, *newv, *newwork;

/* precalc some boundaries */
    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    resid = work + 2 * size;

    scale = 2.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    scale = scale + (2.0 / (gridhy * gridhy * L->get_yside() * L->get_yside()));
    scale = scale + (2.0 / (gridhz * gridhz * L->get_zside() * L->get_zside()));
    scale = step / scale;

//    T.trade_images (f_mat, dimx, dimy, dimz, nb_ids, CENTRAL_TRADE);
    T->trade_images (f_mat, dimx, dimy, dimz, FULL_TRADE);


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
            T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
        }
        else {
            T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
        }

    }


/*
 * on coarsest grid, we are finished
 */

    if (level == max_levels)
    {

        return;

    }                           /* end if */

/* size for next smaller grid */
    dx2 = MG_SIZE (dimx, level, gxsize, gxoffset, pxdim, &ixoff, boundaryflag);
    dy2 = MG_SIZE (dimy, level, gysize, gyoffset, pydim, &iyoff, boundaryflag);
    dz2 = MG_SIZE (dimz, level, gzsize, gzoffset, pzdim, &izoff, boundaryflag);
    siz2 = (dx2 + 2) * (dy2 + 2) * (dz2 + 2);

    // If dx2, dy2 or dz2 is negative then it means that too many multigrid levels were requested so we just return and continue processing.
    // Since this is normally called inside loops we don't print an error message each time but wait until the destructor is called.
    if((dx2 < 0) || (dy2 < 0) || (dz2 < 0)) {
        level_flag++;
        return;
    }

/* evaluate residual */
    eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
    T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);



/* set storage pointers in the current workspace */
    newv = &work[0];
    newf = &work[siz2];
    newwork = &work[2 * siz2];


    for (i = 0; i < mu_cyc; i++)
    {

        mg_restrict (resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

        /* call mgrid solver on new level */
        mgrid_solv(newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                    gridhy * 2.0, gridhz * 2.0, level + 1, nb_ids,
                    max_levels, pre_cyc, post_cyc, mu_cyc, step, k,
                    gxsize, gysize, gzsize,
                    gxoffset, gyoffset, gzoffset,
                    pxdim, pydim, pzdim, boundaryflag);


        mg_prolong (resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        scale = ONE;

        QMD_axpy (size, scale, resid, ione, v_mat, ione);

        /* re-solve on this grid level */

        T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);

        for (cycl = 0; cycl < post_cyc[level]; cycl++)
        {

            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, k);

            /* trade boundary info */
            T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);

        }                       /* end for */

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid);
            T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);

        }                       /* end if */


    }                           /* for mu_cyc */

}

template <typename RmgType>
void Mgrid::mg_restrict (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    int ix, iy, iz, ibrav;
    int incy, incx, incy2, incx2;
    int x0, xp, xm, y0, yp, ym, z0, zp, zm;
    double scale, face, corner, edge;

    ibrav = L->get_ibrav_type();

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    incy2 = dz2 + 2;
    incx2 = (dz2 + 2) * (dy2 + 2);


    switch (ibrav)
    {

        case CUBIC_PRIMITIVE:
        case CUBIC_FC:
        case ORTHORHOMBIC_PRIMITIVE:
        case HEXAGONAL:

            scale = ONE / 64.0;
            for (ix = 1; ix <= dx2; ix++)
            {

                x0 = 2 * ix - 1 + xoffset;
                xp = x0 + 1;
                xm = x0 - 1;

                for (iy = 1; iy <= dy2; iy++)
                {

                    y0 = 2 * iy - 1 + yoffset;
                    yp = y0 + 1;
                    ym = y0 - 1;

                    for (iz = 1; iz <= dz2; iz++)
                    {

                        z0 = 2 * iz - 1 + zoffset;
                        zp = z0 + 1;
                        zm = z0 - 1;

                        face = full[xm * incx + y0 * incy + z0] +
                            full[xp * incx + y0 * incy + z0] +
                            full[x0 * incx + ym * incy + z0] +
                            full[x0 * incx + yp * incy + z0] +
                            full[x0 * incx + y0 * incy + zm] + full[x0 * incx + y0 * incy + zp];

                        corner =
                            full[xm * incx + ym * incy + zm] +
                            full[xm * incx + ym * incy + zp] +
                            full[xm * incx + yp * incy + zm] +
                            full[xm * incx + yp * incy + zp] +
                            full[xp * incx + ym * incy + zm] +
                            full[xp * incx + ym * incy + zp] +
                            full[xp * incx + yp * incy + zm] + full[xp * incx + yp * incy + zp];

                        edge = full[xm * incx + y0 * incy + zm] +
                            full[xm * incx + ym * incy + z0] +
                            full[xm * incx + yp * incy + z0] +
                            full[xm * incx + y0 * incy + zp] +
                            full[x0 * incx + ym * incy + zm] +
                            full[x0 * incx + yp * incy + zm] +
                            full[x0 * incx + ym * incy + zp] +
                            full[x0 * incx + yp * incy + zp] +
                            full[xp * incx + y0 * incy + zm] +
                            full[xp * incx + ym * incy + z0] +
                            full[xp * incx + yp * incy + z0] + full[xp * incx + y0 * incy + zp];


                        half[ix * incx2 + iy * incy2 + iz] =
                            scale * (8.0 * full[x0 * incx + y0 * incy + z0] + 4.0 * face + 2.0 * edge +
                                     corner);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;


        case CUBIC_BC:

            scale = ONE / 52.0;

            for (ix = 1; ix <= dx2; ix++)
            {

                x0 = 2 * ix - 1 + xoffset;
                xp = x0 + 1;
                xm = x0 - 1;

                for (iy = 1; iy <= dy2; iy++)
                {

                    y0 = 2 * iy - 1 + yoffset;
                    yp = y0 + 1;
                    ym = y0 - 1;

                    for (iz = 1; iz <= dz2; iz++)
                    {

                        z0 = 2 * iz - 1 + zoffset;
                        zp = z0 + 1;
                        zm = z0 - 1;

                        face = full[xm * incx + ym * incy + z0] +
                            full[xm * incx + y0 * incy + zm] +
                            full[x0 * incx + ym * incy + zm] +
                            full[x0 * incx + yp * incy + zp] +
                            full[xp * incx + y0 * incy + zp] + full[xp * incx + yp * incy + z0];

                        corner =
                            full[xm * incx + ym * incy + zm] +
                            full[xm * incx + y0 * incy + z0] +
                            full[x0 * incx + ym * incy + z0] +
                            full[x0 * incx + y0 * incy + zm] +
                            full[x0 * incx + y0 * incy + zp] +
                            full[x0 * incx + yp * incy + z0] +
                            full[xp * incx + y0 * incy + z0] + full[xp * incx + yp * incy + zp];


                        half[ix * incx2 + iy * incy2 + iz] =
                            scale * (8.0 * full[x0 * incx + y0 * incy + z0] + 4.0 * corner +
                                     2.0 * face);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        case 20:

            scale = ONE / 80.0;
            for (ix = 1; ix <= dx2; ix++)
            {

                x0 = 2 * ix - 1 + xoffset;
                xp = x0 + 1;
                xm = x0 - 1;

                for (iy = 1; iy <= dy2; iy++)
                {

                    y0 = 2 * iy - 1 + yoffset;
                    yp = y0 + 1;
                    ym = y0 - 1;

                    for (iz = 1; iz <= dz2; iz++)
                    {

                        z0 = 2 * iz - 1 + zoffset;
                        zp = z0 + 1;
                        zm = z0 - 1;

                        face = full[xm * incx + ym * incy + z0] +
                            full[xm * incx + y0 * incy + zm] +
                            full[x0 * incx + ym * incy + zm] +
                            full[x0 * incx + yp * incy + zp] +
                            full[xp * incx + y0 * incy + zp] + full[xp * incx + yp * incy + z0];

                        edge =
                            full[xm * incx + y0 * incy + z0] +
                            full[xm * incx + y0 * incy + zp] +
                            full[xm * incx + yp * incy + z0] +
                            full[x0 * incx + ym * incy + z0] +
                            full[x0 * incx + ym * incy + zp] +
                            full[x0 * incx + y0 * incy + zm] +
                            full[x0 * incx + y0 * incy + zp] +
                            full[x0 * incx + yp * incy + zm] +
                            full[x0 * incx + yp * incy + z0] +
                            full[xp * incx + ym * incy + z0] +
                            full[xp * incx + y0 * incy + zm] + full[xp * incx + y0 * incy + z0];


                        half[ix * incx2 + iy * incy2 + iz] =
                            scale * (8.0 * full[x0 * incx + y0 * incy + z0] + 5.0 * edge + 2.0 * face);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not programmed");

    }                           /* end switch */


}                               /* end mg_restrict */

template <typename RmgType>
void Mgrid::mg_prolong (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    int ix, iy, iz;
    int incx, incy, incxr, incyr;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    incyr = dz2 + 2;
    incxr = (dz2 + 2) * (dy2 + 2);


    /* transfer coarse grid points to fine grid along with the
     * high side image point
     */

    for (ix = 1-xoffset; ix <= dimx/2 + 1; ix++)
    {

        for (iy = 1-yoffset; iy <= dimy/2 + 1; iy++)
        {

            for (iz = 1-zoffset; iz <= dimz/2 + 1; iz++)
            {

                full[(2 * ix - 1+xoffset) * incx + (2 * iy - 1+yoffset) * incy + 2 * iz - 1+zoffset] =
                    half[ix * incxr + iy * incyr + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    /* interior center points
     */
    for (ix = 2-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 2-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 2-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                    0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz + 1] +
                    0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz - 1] +
                    0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz + 1] +
                    0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz - 1] +
                    0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz + 1] +
                    0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz - 1] +
                    0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 1-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 2-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[ix * incx + iy * incy + iz - 1] +
                    0.5 * full[ix * incx + iy * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 2-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 1-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[ix * incx + (iy - 1) * incy + iz] +
                    0.5 * full[ix * incx + (iy + 1) * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 1-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 1-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[(ix - 1) * incx + iy * incy + iz] +
                    0.5 * full[(ix + 1) * incx + iy * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 2-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 2-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[ix * incx + (iy - 1) * incy + iz - 1] +
                    0.25 * full[ix * incx + (iy - 1) * incy + iz + 1] +
                    0.25 * full[ix * incx + (iy + 1) * incy + iz - 1] +
                    0.25 * full[ix * incx + (iy + 1) * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 1-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 2-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[(ix - 1) * incx + iy * incy + iz - 1] +
                    0.25 * full[(ix - 1) * incx + iy * incy + iz + 1] +
                    0.25 * full[(ix + 1) * incx + iy * incy + iz - 1] +
                    0.25 * full[(ix + 1) * incx + iy * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2-xoffset; ix <= dimx; ix += 2)
    {

        for (iy = 2-yoffset; iy <= dimy; iy += 2)
        {

            for (iz = 1-zoffset; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[(ix - 1) * incx + (iy - 1) * incy + iz] +
                    0.25 * full[(ix + 1) * incx + (iy - 1) * incy + iz] +
                    0.25 * full[(ix - 1) * incx + (iy + 1) * incy + iz] +
                    0.25 * full[(ix + 1) * incx + (iy + 1) * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end mg_prolong */



template <typename RmgType>
void Mgrid::eval_residual (RmgType * mat, RmgType * f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, RmgType * res)
{
    int size, idx;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        res[idx] = 0.0;
    FD.app_del2c (mat, res, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    for (idx = 0; idx < size; idx++)
        res[idx] = f_mat[idx] - res[idx];


}                               /* end eval_residual */



template <typename RmgType>
void Mgrid::solv_pois (RmgType * vmat, RmgType * fmat, RmgType * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double k)
{
    int size, idx;
    double scale;
    double diag;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        work[idx] = ZERO;
    diag = -FD.app_del2c (vmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    scale = step / diag;
    
    // Non-zero k effectively means we are solving the Helmholtz rather than Poissons equation
    if(k != 0.0) {

        for (idx = 0; idx < size; idx++)
        {

            vmat[idx] += scale * (work[idx] - k*vmat[idx] - fmat[idx]);

        }                           /* end for */

     }
     else {

        for (idx = 0; idx < size; idx++)
        {

            vmat[idx] += scale * (work[idx] - fmat[idx]);

        }                           /* end for */

     }

}                               /* end solv_pois */



/** Compute 1-D grid sizes for the next multigrid level 

Inputs:<br>
curdim        = current size of this grid on this node<br>
global_dim    = global grid dimension<br>
global_offset = offset of edge of this node grid on the global grid<br>
global_pdim   = dimension of this node grid<br>
bctype        = boundary condition<br>


Outputs:<br>
*roffset      = pointer to grid offset (always 0 or 1)<br>

Return value  = size of next grid level


*/

int Mgrid::MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
{
    int skip, new_dim, istart, istop;

    // Default offset is 0
    *roffset = 0;

    if(bctype == PERIODIC) {

        skip = (2 << curlevel);
        // First check if we have too many multigrid levels. For periodic boundary
        // conditions the next level of the global grid must be divisible by 2
        if ((global_dim % skip) != 0) {
            return -1;
        }

        // Require at least one point in the level
        new_dim = global_pdim / skip;
        if(!new_dim) {
            return -1;
        }

        // evenly divisible then we are done
        if(!(global_pdim % skip)) return new_dim;

        // Check if first point is included and if not subtract
        istart = skip - global_offset % skip;
        istop = (global_offset + global_pdim - 1) % skip;
        if((istart == skip) || (istop == skip)) new_dim++;
        
        // Perform offset check
        if((istart == skip) || (istart == 0)) {
            return new_dim;
        }
        *roffset = 1;
        
        return new_dim;

    }

    rmg_error_handler (__FILE__, __LINE__, "Boundary condition not programmed."); 

}

