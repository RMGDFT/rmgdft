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
#include "packfuncs.h"
#include "boundary_conditions.h"

template void Mgrid::mgrid_solv<float>(float*, float*, float*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<double>(double*, double*, double*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv<float>(float*, float*, float*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int, double *);

template void Mgrid::mgrid_solv<double>(double*, double*, double*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int, double *);

template void Mgrid::mgrid_solv<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int, double *);

template void Mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double, double, double *, int, int, int, int, int, int, int, int, int, int, double *);

template void Mgrid::mgrid_solv_pois<float>(float*, float*, float*, int, int, int, double, double, double, int, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<double>(double*, double*, double*, int, int, int, double, double, double, int, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_pois<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int, int*, int*, int, double, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<float>(float*, float*, float*, int, int, int, double, double, double, int, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<double>(double*, double*, double*, int, int, int, double, double, double, int, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);

template void Mgrid::mgrid_solv_schrodinger<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int, int*, int*, int, double, double *, int, int, int, int, int, int, int, int, int, int);


template void Mgrid::eval_residual (double *, double *, int, int, int, double, double, double, double *, double *);

template void Mgrid::eval_residual (float *, float *, int, int, int, double, double, double, float *, double *);

template void Mgrid::solv_pois (double *, double *, double *, int, int, int, double, double, double, double, double, double, double *);

template void Mgrid::solv_pois (float *, float *, float *, int, int, int, double, double, double, double, double, double, double *);

template void Mgrid::mg_restrict(float*, float*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_restrict(double*, double*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_restrict(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_restrict(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);

template void Mgrid::mg_prolong(float*, float*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong(double*, double*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);

template void Mgrid::mg_prolong_cubic(float*, float*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong_cubic(double*, double*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong_cubic(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void Mgrid::mg_prolong_cubic(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);

//template void Mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, double, double, double, int, int*, int, int*, int*, int, double, double, int, int, int, int, int, int, int, int, int, int);

int Mgrid::level_warning;

Mgrid::Mgrid(Lattice *lptr, TradeImages *tptr)
{
    L = lptr;
    T = tptr;
    level_flag = 0;
    this->ibrav = L->get_ibrav_type();
    this->timer_mode = false;
    this->central_trade = false;
    if((this->ibrav == CUBIC_PRIMITIVE) || (this->ibrav == ORTHORHOMBIC_PRIMITIVE)) this->central_trade = true;
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
                 int level, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, 
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    Mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 gridhx, gridhy, gridhz,
                 level, max_levels, pre_cyc,
                 post_cyc, mu_cyc, step, 0.0, 0.0, NULL,
                 gxsize, gysize, gzsize,
                 gxoffset, gyoffset, gzoffset,
                 pxdim, pydim, pzdim, boundaryflag);

}


template <typename RmgType>
void Mgrid::mgrid_solv_schrodinger (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double *pot,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    Mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 gridhx, gridhy, gridhz,
                 level, max_levels, pre_cyc,
                 post_cyc, mu_cyc, step, 0.0, 0.0, pot,
                 gxsize, gysize, gzsize,
                 gxoffset, gyoffset, gzoffset,
                 pxdim, pydim, pzdim, boundaryflag);

}

template <typename RmgType>
void Mgrid::mgrid_solv (RmgType * __restrict__ v_mat, RmgType * __restrict__ f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double Zfac, double k, double *pot,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag)
{
    double kvec[3] = {0.0, 0.0, 0.0};
    Mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 gridhx, gridhy, gridhz,
                 level, max_levels, pre_cyc,
                 post_cyc, mu_cyc, step, Zfac, k, pot,
                 gxsize, gysize, gzsize,
                 gxoffset, gyoffset, gzoffset,
                 pxdim, pydim, pzdim, boundaryflag, kvec);
}

template <typename RmgType>
void Mgrid::mgrid_solv (RmgType * __restrict__ v_mat, RmgType * __restrict__ f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double Zfac, double k, double *pot,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag, double *kvec)
{
    RmgTimer *RT = NULL;
    BaseThread *Threads = BaseThread::getBaseThread(0);
    RmgType half(0.5);

    std::string timername;
    if(this->timer_mode) {
        timername = "Mgrid_solv: level " + boost::to_string(level);
        RT = new RmgTimer(timername.c_str());
    }

    int ixoff, iyoff, izoff;
    int mindim = std::min(dimx, dimy);
    mindim = std::min(mindim, dimz);

/* precalc some boundaries */
    int size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    RmgType *resid = work + 2 * size;

    double scale = 2.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    scale = scale + (2.0 / (gridhy * gridhy * L->get_yside() * L->get_yside()));
    scale = scale + (2.0 / (gridhz * gridhz * L->get_zside() * L->get_zside()));
    scale = 1.0 / (scale + Zfac);
    scale = step * scale;

    //RmgType *f_mat_t = (RmgType *)alloca(size*sizeof(RmgType));
    RmgType *f_mat_t = (RmgType *)&work[3*size];
    for(int idx=0;idx<size;idx++)f_mat_t[idx] = f_mat[idx];

    bool check = (dimx >= 3) && (dimy >= 3) && (dimz >= 3);
    int tid = Threads->get_thread_tid();
    int active_threads = 1;
    if(tid >= 0) active_threads = Threads->barrier->barrier_count();
    if(active_threads < 2) check = false;
    if(pot || (k != 0.0) || (pre_cyc[level] > MAX_TRADE_IMAGES) || !check)
    {

        T->trade_images (f_mat, dimx, dimy, dimz, FULL_TRADE);
        if(pot) T->trade_images (pot, dimx, dimy, dimz, FULL_TRADE);
//        for (int idx = 0; idx < size; idx++) v_mat[idx] = -(RmgType)scale * f_mat[idx];
        for (int idx = 0; idx < size; idx++) v_mat[idx] = half*(RmgType)scale * f_mat[idx];

        // solve on this grid level 
        for (int cycl = 0; cycl < pre_cyc[level]; cycl++)
        {
            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, k, pot);

            /* trade boundary info */
            if (((level >= max_levels) && (cycl == pre_cyc[level]-1)) || !this->central_trade) {
                T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
            }
            else {
                T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
            }
        }
    }
    else
    {
        // Convert f_mat into p-type work grid then trade images up to 4
        int offset = std::min(mindim, 4);  // offset now holds the max number we can process at once
        CPP_pack_stop (f_mat, work, dimx, dimy, dimz);
        T->trade_imagesx (work, f_mat, dimx, dimy, dimz, offset, FULL_TRADE);
        int fullsteps = pre_cyc[level] / offset;
        int rem = pre_cyc[level] % offset;
        for(int steps=0;steps < fullsteps;steps++)
        {
            size = (dimx + 2*offset)*(dimy + 2*offset)*(dimz + 2*offset);
            if(steps == 0) 
            {
                for (int idx = 0; idx < size; idx++) v_mat[idx] = half*(RmgType)scale * f_mat[idx];
            }
            else
            {
                CPP_pack_stop (v_mat, work, dimx, dimy, dimz);
                T->trade_imagesx (work, v_mat, dimx, dimy, dimz, offset, FULL_TRADE);
            }
            solv_pois_offset (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, offset, offset);
        }
        if(rem)
        {
            CPP_pack_stop (v_mat, work, dimx, dimy, dimz);
            T->trade_imagesx (work, v_mat, dimx, dimy, dimz, rem, FULL_TRADE);
            solv_pois_offset (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, rem, offset);
        }
        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
    }

/*
 * on coarsest grid, we are finished
 */

    if (level >= max_levels)
    {
        if(this->timer_mode) delete RT;
        return;

    }                           /* end if */

/* size for next smaller grid */
    int dx2 = MG_SIZE (dimx, level, gxsize, gxoffset, pxdim, &ixoff, boundaryflag);
    int dy2 = MG_SIZE (dimy, level, gysize, gyoffset, pydim, &iyoff, boundaryflag);
    int dz2 = MG_SIZE (dimz, level, gzsize, gzoffset, pzdim, &izoff, boundaryflag);

    // If dx2, dy2 or dz2 is negative then it means that too many multigrid levels were requested so we just return and continue processing.
    // Since this is normally called inside loops we don't print an error message each time but wait until the destructor is called.
    if((dx2 < 0) || (dy2 < 0) || (dz2 < 0)) {
        level_flag++;
        if(this->timer_mode) delete RT;
        return;
    }


/* evaluate residual */
    eval_residual (v_mat, f_mat_t, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid, pot);
    T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);


/* set storage pointers in the current workspace */
    RmgType *newv = &work[0];
    RmgType *newf = &work[size];
    RmgType *newwork = &work[2 * size];
    double *newpot=NULL;
    if(pot) newpot = &pot[size];


    for (int i = 0; i < mu_cyc; i++)
    {

        mg_restrict (resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        if(pot) mg_restrict (pot, newpot, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

        /* call mgrid solver on new level */
        mgrid_solv(newv, newf, newwork, dx2, dy2, dz2, gridhx * 2.0,
                    gridhy * 2.0, gridhz * 2.0, level + 1,
                    max_levels, pre_cyc, post_cyc, mu_cyc, step, 2.0*Zfac, k, newpot,
                    gxsize, gysize, gzsize,
                    gxoffset, gyoffset, gzoffset,
                    pxdim, pydim, pzdim, boundaryflag, kvec);

        mg_prolong (resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        for(int idx = 0;idx < size;idx++) v_mat[idx] += resid[idx];


        /* re-solve on this grid level */
        if(pot || (k != 0.0) || (pre_cyc[level] > MAX_TRADE_IMAGES) || !check)
        {
            if(!this->central_trade)
                T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
            else
                T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);

            for (int cycl = 0; cycl < post_cyc[level]; cycl++)
            {

                /* solve once */
                solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, k, pot);

                /* trade boundary info */
                if(cycl < (post_cyc[level] - 1))
                {
                    if(!this->central_trade)
                        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
                    else
                        T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
                }

            }                       /* end for */
        }
        else
        {
            int offset = std::min(mindim, 4);  // offset now holds the max number we can process at once
            int fullsteps = post_cyc[level] / offset;
            int rem = post_cyc[level] % offset;
            for(int steps=0;steps < fullsteps;steps++)
            {
                CPP_pack_stop (v_mat, work, dimx, dimy, dimz);
                T->trade_imagesx (work, v_mat, dimx, dimy, dimz, offset, FULL_TRADE);
                solv_pois_offset (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, offset, offset);
            }
            if(rem)
            {
                CPP_pack_stop (v_mat, work, dimx, dimy, dimz);
                T->trade_imagesx (work, v_mat, dimx, dimy, dimz, rem, FULL_TRADE);
                solv_pois_offset (v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, Zfac, rem, offset);
            }
        }
        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);

        /* evaluate max residual */
        if (i < (mu_cyc - 1))
        {

            eval_residual (v_mat, f_mat_t, dimx, dimy, dimz, gridhx, gridhy, gridhz, resid, pot);
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
        case CUBIC_BC:
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
void Mgrid::mg_prolong_cubic (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    RmgType cc[10][4];
    double frac;

    size_t alloc = (dx2 + 4)*(dy2 + 4)*(dz2 + 4);
    RmgType *half_c = new RmgType[alloc];
    RmgType *rptr = new RmgType[alloc];

    int ifxs = (dimy + 2) * (dimz + 2);
    int ifys = dimz + 2;

    int icxs = (dy2 + 4) * (dz2 + 4);
    int icys = dz2 + 4;


    for (int i = 0; i < 2; i++)
    {
        frac = (double) i / (double) 2;
        cc[i][0] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc[i][1] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc[i][2] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc[i][3] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
    }


    //CPP_pack_stop (half, rptr, dx2, dy2, dz2);
    T->trade_imagesx (half, half_c, dx2, dy2, dz2, 2, FULL_TRADE);

    for (int i = 2; i < dx2 + 2; i++)
    {
        for (int j = 2; j < dy2 + 2; j++)
        {
            for (int k = 2; k < dz2 + 2; k++)
            {
                full[(2 * (i - 2) + 1)*ifxs + (2 * (j - 2) + 1)*ifys + 2 * (k - 2) + 1] =
                    half_c[i*icxs + j*icys + k];
            }
        }
    }

    for (int i = 2; i < dx2 + 2; i++)
    {
        for (int j = 2; j < dy2 + 2; j++)
        {
            for (int k = 2; k < dz2 + 2; k++)
            {

                for (int in = 1; in < 2; in++)
                {
                    RmgType tmp1(0.0);
                    RmgType tmp2(0.0);
                    RmgType tmp3(0.0);
                    int basis1 = -1;

                    for (int ii = 0; ii < 4; ii++)
                    {
                        tmp1 += cc[in][ii] * half_c[(i + basis1)*icxs + j*icys + k];
                        tmp2 += cc[in][ii] * half_c[i*icxs + (j + basis1)*icys + k];
                        tmp3 += cc[in][ii] * half_c[i*icxs + j*icys + k + basis1];
                        ++basis1;
                    }

                    full[(2 * (i - 2) + in + 1)*ifxs + (2 * (j - 2) + 1)*ifys + 2 * (k - 2) + 1] = tmp1;
                    full[(2 * (i - 2) + 1)*ifxs + (2 * (j - 2) + in + 1)*ifys + 2 * (k - 2) + 1] = tmp2;
                    full[(2 * (i - 2) + 1)*ifxs + (2 * (j - 2) + 1)*ifys + 2 * (k - 2) + in + 1] = tmp3;
                }

            }
        }
    }

    for (int i = 2; i < dx2 + 2; i++)
    {
        for (int j = 2; j < dy2 + 2; j++)
        {
            for (int k = 2; k < dz2 + 2; k++)
            {

                for (int in = 1; in < 2; in++)
                {
                    for (int jn = 1; jn < 2; jn++)
                    {

                        RmgType tmp1(0.0);
                        RmgType tmp2(0.0);
                        RmgType tmp3(0.0);
                        int basis1 = -1;

                        for (int ii = 0; ii < 4; ii++)
                        {
                            int basis2 = -1;
                            for (int jj = 0; jj < 4; jj++)
                            {
                                tmp1 +=
                                    cc[in][ii] * cc[jn][jj] * half_c[(i + basis1)*icxs +(j + basis2)*icys + k];
                                tmp2 +=
                                    cc[in][ii] * cc[jn][jj] * half_c[(i + basis1)*icxs + j*icys + k + basis2];
                                tmp3 +=
                                    cc[in][ii] * cc[jn][jj] * half_c[i*icxs + (j + basis1)*icys + k + basis2];
                                ++basis2;
                            }
                            ++basis1;
                        }
                        full[(2 * (i - 2) + in + 1)*ifxs + (2 * (j - 2) + jn + 1)*ifys + 2 * (k - 2) + 1] = tmp1;
                        full[(2 * (i - 2) + in + 1)*ifxs + (2 * (j - 2) + 1)*ifys + 2 * (k - 2) + jn + 1] = tmp2;
                        full[(2 * (i - 2) + 1)*ifxs + (2 * (j - 2) + in + 1)*ifys + 2 * (k - 2) + jn + 1] = tmp3;
                    }
                }
            }
        }
    }


    for (int i = 2; i < dx2 + 2; i++)
    {
        for (int j = 2; j < dy2 + 2; j++)
        {
            for (int k = 2; k < dz2 + 2; k++)
            {

                for (int in = 1; in < 2; in++)
                {
                    for (int jn = 1; jn < 2; jn++)
                    {
                        for (int kn = 1; kn < 2; kn++)
                        {

                            RmgType tmp1(0.0);
                            int basis1 = -1;
                            for (int ii = 0; ii < 4; ii++)
                            {
                                int basis2 = -1;
                                for (int jj = 0; jj < 4; jj++)
                                {
                                    int basis3 = -1;
                                    for (int kk = 0; kk < 4; kk++)
                                    {
                                        tmp1 +=
                                            cc[in][ii] * cc[jn][jj] * cc[kn][kk] * half_c[(i + basis1)*icxs + (j + basis2)*icys + k + basis3];
                                        ++basis3;
                                    }
                                    ++basis2;
                                }
                                ++basis1;
                            }
                            full[(2 * (i - 2) + in + 1)*ifxs + (2 * (j - 2) + jn + 1)*ifys + 2 * (k - 2) + kn + 1] = tmp1;

                        }
                    }
                }
            }
        }
    }

    delete [] rptr;
    delete [] half_c;
}

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
//    for (idx = 0; idx < size; idx++) work[idx] = 0.0;
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


// Used to handle multiple sweeps case. By using a higher level trade images once latency is reduced
// at the cost of doing more local work.
//
// vmat is shrunk on each step while fmat is kept fixed
template <typename RmgType>
void Mgrid::solv_pois_offset (RmgType * __restrict__ vmat, RmgType * __restrict__ fmat, RmgType * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double Zfac, int offset, int foffset)
{
    FiniteDiff FD(L);
    int incy = (dimz + 2*foffset);
    int incx = (dimy + 2*foffset)*(dimz + 2*foffset);

    while(offset > 0)
    {
        offset--;
        double diag = -FD.app2_del2 (vmat, work, dimx+2*offset, dimy+2*offset, dimz+2*offset, gridhx, gridhy, gridhz);
        double scale = 1.0 / (diag + Zfac);
        scale = step * scale;
        int idx = 0;
        int sidx = foffset - offset - 1;
        int eidx = foffset + offset + 1;
        for(int ix=sidx;ix<eidx+dimx;ix++)
        {
            for(int iy=sidx;iy<eidx+dimy;iy++)
            {
                for(int iz=sidx;iz<eidx+dimz;iz++)
                {
                    vmat[idx] += (RmgType)scale * (work[idx] - fmat[ix*incx + iy*incy + iz]);
                    idx++;
                }
            }
        }
        if(offset)
        {
            CPP_pack_stop (vmat, vmat, dimx+2*offset, dimy+2*offset, dimz+2*offset);
        }
    }
}                               /* end solv_pois_offset */


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

