/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/



#include <complex>
#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "Mgrid.h"
#include "FiniteDiff.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "boundary_conditions.h"


template void Mgrid::mgrid_solv<float>(float*, float*, float*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<double>(double*, double*, double*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<float>(float*, float*, float*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<double>(double*, double*, double*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<float>(float*, float*, float*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<double>(double*, double*, double*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);


template void Mgrid::eval_residual (double *, double *, int, int, int, double, double, double, double *, double *);

template void Mgrid::eval_residual (float *, float *, int, int, int, double, double, double, float *, double *);

template void Mgrid::solv_pois (double *, double *, double *, int, int, int, double, double, double, double, double, double, double *);

template void Mgrid::solv_pois (float *, float *, float *, int, int, int, double, double, double, double, double, double, double *);

//template void Mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, int, int, int, int, int, int, int, int, int, int);

int Mgrid::level_warning;

Mgrid::Mgrid(Lattice *lptr, TradeImages *tptr)
{
    L = lptr;
    T = tptr;
    level_flag = 0;
    this->timer_mode = false;
}

Mgrid::~Mgrid(void)
{
    if(level_flag && !Mgrid::level_warning)
        std::cout << "Warning: too many multigrid levels were requested " << level_flag << " times.\n";
    Mgrid::level_warning = true;   // Only want to print one warning
}

void Mgrid::set_timer_mode(bool verbose)
{
    Mgrid::timer_mode = verbose;
}


// Poisson variant
template <typename RmgType>
void Mgrid::mgrid_solv_pois (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, 
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    Mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 gridhx, gridhy, gridhz,
                 level, nb_ids, max_levels, pre_cyc,
                 post_cyc, mu_cyc, step, 0.0, 0.0, NULL,
                 gxsize, gysize, gzsize,
                 gxoffset, gyoffset, gzoffset,
                 pxdim, pydim, pzdim, boundaryflag);

}


template <typename RmgType>
void Mgrid::mgrid_solv_schrodinger (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double *pot,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    Mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 gridhx, gridhy, gridhz,
                 level, nb_ids, max_levels, pre_cyc,
                 post_cyc, mu_cyc, step, 0.0, 0.0, pot,
                 gxsize, gysize, gzsize,
                 gxoffset, gyoffset, gzoffset,
                 pxdim, pydim, pzdim, boundaryflag);

}


template <typename RmgType>
void Mgrid::mgrid_solv (RmgType * __restrict__ v_mat, RmgType * __restrict__ f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double Zfac, double k, double *pot,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    RmgTimer *RT = NULL;
    std::string timername;
    if(this->timer_mode) {
        timername = "Mgrid_solv: level " + boost::to_string(level);
        RT = new RmgTimer(timername.c_str());
    }

    int i;
    int cycl;
    int size, idx;
    double scale;
    int dx2, dy2, dz2, siz2;
    int ixoff, iyoff, izoff;
    RmgType *resid, *newf, *newv, *newwork;
    double *newpot=NULL;

/* precalc some boundaries */
    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    resid = work + 2 * size;

    scale = 2.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    scale = scale + (2.0 / (gridhy * gridhy * L->get_yside() * L->get_yside()));
    scale = scale + (2.0 / (gridhz * gridhz * L->get_zside() * L->get_zside()));
    scale = 1.0 / (scale + Zfac);
    scale = step * scale;


//    T.trade_images (f_mat, dimx, dimy, dimz, nb_ids, CENTRAL_TRADE);
    T->trade_images (f_mat, dimx, dimy, dimz, FULL_TRADE);


    if(pot) 
    {
        T->trade_images (pot, dimx, dimy, dimz, FULL_TRADE);
        for (idx = 0; idx < size; idx++) v_mat[idx] = -(RmgType)scale * (f_mat[idx] + (RmgType)pot[idx]*v_mat[idx]);
    }
    for (idx = 0; idx < size; idx++) v_mat[idx] = -(RmgType)scale * f_mat[idx];


/*
 * solve on this grid level 
 */


    for (cycl = 0; cycl < pre_cyc[level]; cycl++)
    {

        /* solve once */
        solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, k, pot);


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
        if(this->timer_mode) delete RT;
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
        if(this->timer_mode) delete RT;
        return;
    }

/* evaluate residual */
    eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid, pot);
    T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);



/* set storage pointers in the current workspace */
    newv = &work[0];
    newf = &work[siz2];
    if(pot) {
        newpot = &pot[siz2];
    }
    else {
        newpot = NULL;
    }
    newwork = &work[2 * siz2];


    for (i = 0; i < mu_cyc; i++)
    {

        mg_restrict (resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        if(pot) mg_restrict (pot, newpot, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

        /* call mgrid solver on new level */
        mgrid_solv(newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                    gridhy * 2.0, gridhz * 2.0, level + 1, nb_ids,
                    max_levels, pre_cyc, post_cyc, mu_cyc, step, 2.0*Zfac, k, newpot,
                    gxsize, gysize, gzsize,
                    gxoffset, gyoffset, gzoffset,
                    pxdim, pydim, pzdim, boundaryflag);


        mg_prolong (resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        scale = 1.0;

        for(int idx = 0;idx < size;idx++)
        {
            v_mat[idx] += resid[idx];
        }


        /* re-solve on this grid level */

        T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);

        for (cycl = 0; cycl < post_cyc[level]; cycl++)
        {

            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, k, pot);

            /* trade boundary info */
            if(cycl < (post_cyc[level] - 1))
                T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
            else
                T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);

        }                       /* end for */

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            eval_residual (v_mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid, pot);
            T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);

        }                       /* end if */

    }                           /* for mu_cyc */

    if(this->timer_mode) delete RT;
}

template <typename RmgType>
void Mgrid::mg_restrict (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    int ix, iy, iz, ibrav;
    int incy, incx, incy2, incx2;
    int x0, xp, xm, y0, yp, ym, z0, zp, zm;
    double scale;
    RmgType face, corner, edge;
    RmgType inplane, zaxis, outplane;

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

            scale = 1.0 / 64.0;
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
                            (RmgType)scale * ((RmgType)8.0 * full[x0 * incx + y0 * incy + z0] + (RmgType)4.0 * face + (RmgType)2.0 * edge +
                                     corner);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        case HEXAGONAL:

            scale = 1.0 / 32.0;

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

                        inplane =
                            full[x0 * incx + yp * incy + z0] +
                            full[x0 * incx + ym * incy + z0] +
                            full[xm * incx + y0 * incy + z0] +
                            full[xp * incx + y0 * incy + z0] +
                            full[xm * incx + yp * incy + z0] +
                            full[xp * incx + ym * incy + z0];

                        zaxis =
                            full[x0 * incx + y0 * incy + zp] +
                            full[x0 * incx + y0 * incy + zm];

                        outplane =
                            full[x0 * incx + yp * incy + zm] +
                            full[x0 * incx + ym * incy + zm] +
                            full[xm * incx + y0 * incy + zm] +
                            full[xp * incx + y0 * incy + zm] +
                            full[xm * incx + yp * incy + zm] +
                            full[xp * incx + ym * incy + zm] +

                            full[x0 * incx + yp * incy + zp] +
                            full[x0 * incx + ym * incy + zp] +
                            full[xm * incx + y0 * incy + zp] +
                            full[xp * incx + y0 * incy + zp] +
                            full[xm * incx + yp * incy + zp] +
                            full[xp * incx + ym * incy + zp];

                        half[ix * incx2 + iy * incy2 + iz] =
                            (RmgType)scale * ((RmgType)4.0 * full[x0 * incx + y0 * incy + z0] + 
                            (RmgType)2.0 * inplane + (RmgType)2.0 * zaxis + outplane);

                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        case CUBIC_BC:

            scale = 1.0 / 52.0;

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
                            (RmgType)scale * ((RmgType)8.0 * full[x0 * incx + y0 * incy + z0] + (RmgType)4.0 * corner +
                                     (RmgType)2.0 * face);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        case 20:

            scale = 1.0 / 80.0;
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
                            (RmgType)scale * ((RmgType)8.0 * full[x0 * incx + y0 * incy + z0] + (RmgType)5.0 * edge + (RmgType)2.0 * face);


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not programmed");

    }                           /* end switch */


}                               /* end mg_restrict */




template <typename RmgType>
void Mgrid::mg_prolong (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    int ix, iy, iz;
    int incx, incy, incxr, incyr;

    int ibrav = L->get_ibrav_type();

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    incyr = dz2 + 2;
    incxr = (dz2 + 2) * (dy2 + 2);

    switch (ibrav)
    {

        case CUBIC_PRIMITIVE:
        case CUBIC_FC:
        case ORTHORHOMBIC_PRIMITIVE:

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
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz + 1];

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
                            (RmgType)0.5 * full[ix * incx + iy * incy + iz - 1] +
                            (RmgType)0.5 * full[ix * incx + iy * incy + iz + 1];

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
                            (RmgType)0.5 * full[ix * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.5 * full[ix * incx + (iy + 1) * incy + iz];

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
                            (RmgType)0.5 * full[(ix - 1) * incx + iy * incy + iz] +
                            (RmgType)0.5 * full[(ix + 1) * incx + iy * incy + iz];

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
                            (RmgType)0.25 * full[ix * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.25 * full[ix * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.25 * full[ix * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.25 * full[ix * incx + (iy + 1) * incy + iz + 1];

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
                            (RmgType)0.25 * full[(ix - 1) * incx + iy * incy + iz - 1] +
                            (RmgType)0.25 * full[(ix - 1) * incx + iy * incy + iz + 1] +
                            (RmgType)0.25 * full[(ix + 1) * incx + iy * incy + iz - 1] +
                            (RmgType)0.25 * full[(ix + 1) * incx + iy * incy + iz + 1];

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
                            (RmgType)0.25 * full[(ix - 1) * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix + 1) * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix - 1) * incx + (iy + 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix + 1) * incx + (iy + 1) * incy + iz];

                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */

            break;

        case HEXAGONAL:

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
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz + 1];

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
                            (RmgType)0.5 * full[ix * incx + iy * incy + iz - 1] +
                            (RmgType)0.5 * full[ix * incx + iy * incy + iz + 1];

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
                            (RmgType)0.5 * full[ix * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.5 * full[ix * incx + (iy + 1) * incy + iz];

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
                            (RmgType)0.5 * full[(ix - 1) * incx + iy * incy + iz] +
                            (RmgType)0.5 * full[(ix + 1) * incx + iy * incy + iz];

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
                            (RmgType)0.25 * full[ix * incx + (iy - 1) * incy + iz - 1] +
                            (RmgType)0.25 * full[ix * incx + (iy - 1) * incy + iz + 1] +
                            (RmgType)0.25 * full[ix * incx + (iy + 1) * incy + iz - 1] +
                            (RmgType)0.25 * full[ix * incx + (iy + 1) * incy + iz + 1];

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
                            (RmgType)0.25 * full[(ix - 1) * incx + iy * incy + iz - 1] +
                            (RmgType)0.25 * full[(ix - 1) * incx + iy * incy + iz + 1] +
                            (RmgType)0.25 * full[(ix + 1) * incx + iy * incy + iz - 1] +
                            (RmgType)0.25 * full[(ix + 1) * incx + iy * incy + iz + 1];

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
                            (RmgType)0.25 * full[(ix - 1) * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix + 1) * incx + (iy - 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix - 1) * incx + (iy + 1) * incy + iz] +
                            (RmgType)0.25 * full[(ix + 1) * incx + (iy + 1) * incy + iz];

                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */

            break;
        
        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not programmed");

    } // end switch

}                               /* end mg_prolong */



template <typename RmgType>
void Mgrid::eval_residual (RmgType * __restrict__ mat, RmgType * __restrict__ f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, RmgType * res, double *pot)
{
    int size, idx;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        res[idx] = 0.0;
    FD.app2_del2 (mat, res, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    if(pot) {
        for (idx = 0; idx < size; idx++) res[idx] = f_mat[idx] + (RmgType)pot[idx]*mat[idx] - res[idx];
    }
    else {
        for (idx = 0; idx < size; idx++) res[idx] = f_mat[idx] - res[idx];
    }


}                               /* end eval_residual */



template <typename RmgType>
void Mgrid::solv_pois (RmgType * __restrict__ vmat, RmgType * __restrict__ fmat, RmgType * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double Zfac, double k, double *pot)
{
    int size, idx;
    double scale;
    double diag;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (idx = 0; idx < size; idx++)
        work[idx] = 0.0;
    diag = -FD.app2_del2 (vmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz);

    scale = 1.0 / (diag + Zfac);
    scale = step * scale;
 
    // Non-zero k effectively means we are solving the Helmholtz rather than Poissons equation
    if(k != 0.0) {

        for (idx = 0; idx < size; idx++)
        {

            vmat[idx] += (RmgType)scale * (work[idx] - (RmgType)k*vmat[idx] - fmat[idx]);

        }                           /* end for */

     }
     else {

        if(pot) 
        {
            for (idx = 0; idx < size; idx++) vmat[idx] += (RmgType)scale * (work[idx] - vmat[idx]*(RmgType)pot[idx] - fmat[idx]);
        }
        else 
        {
            for (idx = 0; idx < size; idx++) vmat[idx] += (RmgType)scale * (work[idx] - fmat[idx]);
        }

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
    return -1;

}

