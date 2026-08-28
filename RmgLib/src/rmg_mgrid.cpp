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
#include <iostream>
#include <vector>
#include <cmath>
#include "rmg_mgrid.h"
#include "FiniteDiff.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "packfuncs.h"
#include "boundary_conditions.h"
#include "rmg_sum_all.h"

namespace rmg
{

template void mgrid::anchor_residual<double> (int , int, int, double *);
template void mgrid::anchor_residual<float> (int , int, int, float *);
template void mgrid::anchor_residual<std::complex<double>> (int , int, int, std::complex<double> *);
template void mgrid::anchor_residual<std::complex<float>> (int , int, int, std::complex<float> *);
template void mgrid::anchor_residual<double> (int, int , int, int, double *);
template void mgrid::anchor_residual<float> (int, int , int, int, float *);
template void mgrid::anchor_residual<std::complex<double>> (int, int , int, int, std::complex<double> *);
template void mgrid::anchor_residual<std::complex<float>> (int, int , int, int, std::complex<float> *);

std::vector<int> mgrid::toffsets;

// Requires r to be in a p type grid with no ghost points
template <typename RmgType>
void mgrid::anchor_residual(int id, int level, int n, RmgType *r)
{

    double s1[2]{0.0,0.0};
    for(int i=0;i < n;i++)
    {
        s1[0] += std::real(r[i]);
        s1[1] += std::imag(r[i]);
    }
    s1[0] = rmg::sum_all(s1[0], this->T->comm);
    s1[1] = rmg::sum_all(s1[1], this->T->comm);
    RmgType scale;
    int global_n = this->gxsize*this->gysize*this->gzsize / std::pow(8, level);
    if constexpr(std::is_same_v<RmgType, std::complex<double>> || std::is_same_v<RmgType, std::complex<float>>)
    {
        RmgType I_t(0.0, 1.0);
        scale = s1[0] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= scale;
        scale = s1[1] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= I_t*scale;
//printf("SSSS  %d  %d  %d  %d  %14.8e  %14.8e\n",dddddd, level, n, global_n, std::real(scale), std::imag(scale));
    }
    else
    {
        scale = s1[0] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= scale;
//printf("SSSS  %d  %d  %d  %14.8e\n", level, n, global_n, scale);
    }
}

// r is in an s type grid with ghost points
template <typename RmgType> void mgrid::anchor_residual(int level, int dimx, int dimy, int dimz, RmgType *r)
{
    double s1[2]{0.0,0.0};

    int incys = dimz + 2;
    int incxs = (dimy + 2) * (dimz + 2);

    for(int ix=1;ix <= dimx;ix++)
    {
        for(int iy=1;iy <= dimy;iy++)
        {
            for(int iz=1;iz <= dimz;iz++)
            {
                s1[0] += std::real(r[ix*incxs + iy*incys + iz]);
                s1[1] += std::imag(r[ix*incxs + iy*incys + iz]);
            }
        }
    }

    // Now we have to scale all the data including the ghost points
    int n = (dimx+2)*(dimy+2)*(dimz+2);
    s1[0] = rmg::sum_all(s1[0], this->T->comm);
    s1[1] = rmg::sum_all(s1[1], this->T->comm);
    RmgType scale;
    int global_n = this->gxsize*this->gysize*this->gzsize / std::pow(8, level);
    if constexpr(std::is_same_v<RmgType, std::complex<double>> || std::is_same_v<RmgType, std::complex<float>>)
    {
        RmgType I_t(0.0, 1.0);
        scale = s1[0] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= scale;
        scale = s1[1] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= I_t*scale;
    }
    else
    {
        scale = s1[0] / (double)global_n;
        for(int i=0;i < n;i++) r[i] -= scale;
    }

    if(1)
    {
        s1[0] = 0.0;
        s1[1] = 0.0;
        for(int ix=1;ix <= dimx;ix++)
        {
            for(int iy=1;iy <= dimy;iy++)
            {
                for(int iz=1;iz <= dimz;iz++)
                {
                    s1[0] += std::real(r[ix*incxs + iy*incys + iz]);
                    s1[1] += std::imag(r[ix*incxs + iy*incys + iz]);
                }
            }
        }

        // Now we have to scale all the data including the ghost points
        s1[0] = rmg::sum_all(s1[0], this->T->comm);
        s1[1] = rmg::sum_all(s1[1], this->T->comm);

    }
}

double mgrid::rel_sradius(int level)
{
    const double PI = 3.14159265358979323;
    int nmax = std::max(std::max(gxsize, gysize), gzsize);
    nmax = nmax / (level + 1);
    return (2.0 + cos(2.0*PI/(double)nmax)) / 3.0;
}

#if 1
//  Computes multigrid smoothing parameters based on Cheybshev polynomials
std::vector<double> pois_periodic_coeffs(
       int nx, int ny, int nz,           // global grid points at this level
       double dx, double dy, double dz,  // grid spacing at this level
       double sigma,                     // smoothing window (0.0 = full band, 1.0 = higher)
       int nsteps)                       // number of steps
{
    const double PI = 3.14159265358979323;
    std::vector<double> coefs(nsteps, 0.0);

    // 1. Directional weights
    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);
    double inv_dz2 = 1.0 / (dz * dz);

    // Gvector units
    double cx = std::cos((2.0 * PI) / nx);
    double cy = std::cos((2.0 * PI) / ny);
    double cz = std::cos((2.0 * PI) / nz);

    double A = (inv_dx2 * cx) + (inv_dy2 * cy) + (inv_dz2 * cz);
    double B = inv_dx2 + inv_dy2 + inv_dz2;
    double r = A / B;

    if (sigma < 0.001 && nx > 4 && ny > 4 && nz > 4) {
        r = 1.0;
    }

    // smoothing window
    double rs = ((1.0 - sigma) * r) / (2.0 + (1.0 + sigma) * r);

    // Chebyshev recursion
    coefs[0] = 1.0;
    for (int k = 1; k < nsteps; ++k) {
        coefs[k] = 4.0 / (4.0 - rs * rs * coefs[k - 1]);
    }

    return coefs;
}
#else
//  Computes multigrid smoothing parameters based on Cheybshev polynomials
std::vector<double> pois_periodic_coeffs(
       int nx, int ny, int nz,           // global grid points at this level
       double dx, double dy, double dz,  // grid spacing at this level
       double sigma_in,                     // smoothing window (0.0 = full band, 1.0 = higher)
       int nsteps)                       // number of steps
{
    std::vector<double> coefs(nsteps, 0.0);
double lmax=12.5;
double lmin=0.25*lmax;
if(sigma_in < 0.5)
{
lmax=12.5;
lmin=0.05;
}
    const double theta = 0.5*(lmax + lmin);
    const double delta = 0.5*(lmax - lmin);
    const double sigma = theta / delta;

    double a=0.0, b=0;
    for(int k=0;k < nsteps;k++)
    {
        if (k==0)
        {
            a = 2.0 / theta; b = 0.0;
        }
        else if (k==1)
        {
            a = 2.0 / theta;
            b  = (a * delta * delta) / 4.0;
        }
        else
        {
            const double beta_new = 1.0 / (2.0*sigma - b);
            b = beta_new;
            a = 2.0 * beta_new / delta;
        }
        coefs[k] = a;
#if 0
        if(is_jacobi)
        {
            a = 2.0/3.0;
            b = 0.0;
        }
#endif
    }
    return coefs;
}
#endif

template void mgrid::mgrid_solv<float>(float*, float*, float*, int, int, int, int, int, double, double, double *, int, int, int);

template void mgrid::mgrid_solv<double>(double*, double*, double*, int, int, int, int, int, double, double, double *, int, int, int);

template void mgrid::mgrid_solv<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, int, int, double, double, double *, int, int, int);

template void mgrid::mgrid_solv<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, int, int, double, double, double *, int, int, int);

template void mgrid::mgrid_solv_pois<float>(float*, float*, float*, int, int, int, int, int, double, int, int, int);

template void mgrid::mgrid_solv_pois<double>(double*, double*, double*, int, int, int, int, int, double, int, int, int);

template void mgrid::mgrid_solv_pois<std::complex <double> >(std::complex<double>*, std::complex<double>*, std::complex<double>*, int, int, int, int, int, double, int, int, int);

template void mgrid::mgrid_solv_pois<std::complex <float> >(std::complex<float>*, std::complex<float>*, std::complex<float>*, int, int, int, int, int, double, int, int, int);


template void mgrid::eval_residual (double *, double *, double *, int, int, int, double, double, double, double *, double *);

template void mgrid::eval_residual (float *, float *, float *, int, int, int, double, double, double, float *, double *);

template void mgrid::solv_pois (double *, double *, double *, int, int, int, double, double, double, double, double, double *);

template void mgrid::solv_pois (float *, float *, float *, int, int, int, double, double, double, double, double, double *);

template void mgrid::mg_restrict(float*, float*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_restrict(double*, double*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_restrict(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_restrict(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);

template void mgrid::mg_prolong(float*, float*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong(double*, double*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);

template void mgrid::mg_prolong_cubic(float*, float*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong_cubic(double*, double*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong_cubic(std::complex<float>*, std::complex<float>*, int, int, int, int, int, int, int, int, int);
template void mgrid::mg_prolong_cubic(std::complex<double>*, std::complex<double>*, int, int, int, int, int, int, int, int, int);


int mgrid::level_warning;

mgrid::mgrid(Lattice *lptr, TradeImages *tptr, rmg::grid *gptr, int density_in, double zmax_in)
{
    L = lptr;
    T = tptr;
    G = gptr;
    density = density_in;
    zmax = zmax_in;

    // Global grid sizes and offsets at the finest level
    gxsize = G->get_NX_GRID(density);
    gysize = G->get_NY_GRID(density);
    gzsize = G->get_NZ_GRID(density);
    gxoffset = G->get_PX_OFFSET(density);
    gyoffset = G->get_PY_OFFSET(density);
    gzoffset = G->get_PZ_OFFSET(density);

    // Grid spacings at all levels
    double hx0 = G->get_hxgrid(density);
    double hy0 = G->get_hygrid(density);
    double hz0 = G->get_hzgrid(density);

    // set up all grid spacings
    int l = 1;
    for (int i=0;i < MAX_MG_LEVELS;i++)
    {
        hx[i] = hx0 * (double)l; 
        hy[i] = hy0 * (double)l; 
        hz[i] = hz0 * (double)l; 
        l *= 2;
    }

    level_flag = 0;
    this->ibrav = L->get_ibrav_type();
    this->timer_mode = false;
    this->central_trade = false;
    if((this->ibrav == CUBIC_PRIMITIVE) || 
       (this->ibrav == ORTHORHOMBIC_PRIMITIVE) || 
       (this->ibrav == TETRAGONAL_PRIMITIVE)) this->central_trade = true;
}

mgrid::~mgrid(void)
{
    if(level_flag && !mgrid::level_warning)
        std::cout << "Warning: too many multigrid levels were requested " << level_flag << " times.\n";
    mgrid::level_warning = true;   // Only want to print one warning
}

void mgrid::set_timer_mode(bool verbose)
{
    mgrid::timer_mode = verbose;
}


// Poisson variant
template <typename RmgType>
void mgrid::mgrid_solv_pois (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 int level, int max_levels, double step, 
                 int pxdim, int pydim, int pzdim)
{
    mgrid::mgrid_solv (v_mat, f_mat, work,
                 dimx, dimy, dimz,
                 level, max_levels, step, 0.0, NULL,
                 pxdim, pydim, pzdim);

}


template <typename RmgType>
void mgrid::mgrid_solv (RmgType * __restrict__ v_mat, RmgType * __restrict__ f_mat, RmgType * work,
                 int dimx, int dimy, int dimz,
                 int level, int max_levels, double step, double k, double *pot,
                 int pxdim, int pydim, int pzdim)
{
    RmgTimer *RT = NULL;
    std::string timername;
    if(this->timer_mode) {
        timername = "mgrid_solv: level " + std::to_string(level);
        RT = new RmgTimer(timername.c_str());
    }

//printf("SSSS  %d   %f\n",level,rel_sradius(level));
    int ixoff, iyoff, izoff;

/* get sizes for the next level */
    int dx2 = MG_SIZE (dimx, level, gxsize, gxoffset, pxdim, &ixoff, boundary_flag);
    int dy2 = MG_SIZE (dimy, level, gysize, gyoffset, pydim, &iyoff, boundary_flag);
    int dz2 = MG_SIZE (dimz, level, gzsize, gzoffset, pzdim, &izoff, boundary_flag);
//printf("LLLLLLLLLL  %d  %d  %d  %d  %d  %d  %d\n",level,dimx,dimy,dimz, pxdim, pydim, pzdim);
    int presweeps = pre_cyc[level];
    std::vector<double> pcoefs;
    int fac = std::pow(2, level);
    int nx = this->T->G->get_NX_GRID(1)/fac;
    int ny = this->T->G->get_NY_GRID(1)/fac;
    int nz = this->T->G->get_NZ_GRID(1)/fac;
    // More sweeps on coarsest level
    if((level >= max_levels) || (dx2 < 0) || (dy2 < 0) || (dz2 < 0))
    {
        int minpts = std::min(std::min(nx, ny), nz);
        int maxpts = std::max(std::max(nx, ny), nz);
        presweeps = std::max(maxpts, 12);
        if(presweeps > minpts) presweeps = minpts;
        pcoefs = pois_periodic_coeffs(nx, ny, nz, hx[level], hy[level], hz[level], 0.0, presweeps);
    }
    else
    {
        pcoefs = pois_periodic_coeffs(nx, ny, nz, hx[level], hy[level], hz[level], 1.0, presweeps);
    }

/* precalc some boundaries */
    int size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    RmgType *resid = work + 2 * size;
    double pscale = std::pow(0.5, level);
    T->trade_images (f_mat, dimx, dimy, dimz, FULL_TRADE);
    if(pot) T->trade_images (pot, dimx, dimy, dimz, FULL_TRADE);
    for (int idx = 0; idx < size; idx++) v_mat[idx] = 0.0;

    // solve on this grid level 
    for (int cycl = 0; cycl < presweeps; cycl++)
    {
        /* solve once */
        solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, hx[level], hy[level], hz[level], pscale*pcoefs[cycl], k, pot);

        /* trade boundary info */
        if (((level >= max_levels) && (cycl == presweeps-1)) || !this->central_trade) {
            T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
        }
        else {
            T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
        }
    }

/*
 * on coarsest grid, we are finished
 */

    if (level >= max_levels)
    {
        anchor_residual(level, dimx, dimy, dimz, v_mat);
        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
        if(this->timer_mode) delete RT;
        return;

    }                           /* end if */

    // If dx2, dy2 or dz2 is negative then it means that too many multigrid levels were requested so we just return and continue processing.
    // Since this is normally called inside loops we don't print an error message each time but wait until the destructor is called.
    if((dx2 < 0) || (dy2 < 0) || (dz2 < 0)) {
        anchor_residual(level, dimx, dimy, dimz, v_mat);
        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
        level_flag++;
        if(this->timer_mode) delete RT;
        return;
    }



/* set storage pointers in the current workspace */
    RmgType *newv = &work[0];
    RmgType *newf = &work[size];
    RmgType *newwork = &work[2 * size];
    double *newpot=NULL;
    if(pot) newpot = &pot[size];


    for (int i = 0; i < mu_cyc[level]; i++)
    {

        /* evaluate residual */
        eval_residual (v_mat, f_mat, work, dimx, dimy, dimz, hx[level], hy[level], hz[level], resid, pot);
        anchor_residual(level, dimx, dimy, dimz, resid);
        T->trade_images (resid, dimx, dimy, dimz, FULL_TRADE);
        mg_restrict (resid, newf, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        if(pot) mg_restrict (pot, newpot, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);

        /* call mgrid solver on new level */
        mgrid_solv(newv, newf, newwork, dx2, dy2, dz2, level + 1,
                    max_levels, step, k, newpot,
                    pxdim, pydim, pzdim);

        mg_prolong (resid, newv, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        for(int idx = 0;idx < size;idx++) v_mat[idx] += resid[idx];

        /* re-solve on this grid level */
        if(!this->central_trade)
            T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
        else
            T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);

        for (int cycl = 0; cycl < post_cyc[level]; cycl++)
        {

            /* solve once */
            solv_pois (v_mat, f_mat, work, dimx, dimy, dimz, hx[level], hy[level], hz[level], pscale*pcoefs[cycl], k, pot);

            /* trade boundary info */
            if(cycl < (post_cyc[level] - 1))
            {
                if(!this->central_trade)
                    T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);
                else
                    T->trade_images (v_mat, dimx, dimy, dimz, CENTRAL_TRADE);
            }

        }                       /* end for */
        anchor_residual(level, dimx, dimy, dimz, v_mat);
        T->trade_images (v_mat, dimx, dimy, dimz, FULL_TRADE);

    }                           /* for mu_cyc */

    if(this->timer_mode) delete RT;
}

template <typename RmgType>
void mgrid::mg_restrict (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
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
        case MONOCLINIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:
        case TRICLINIC_PRIMITIVE:
        case No_Lattice:

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

        case HEXAGONAL2:
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

        case -CUBIC_BC:
        case CUBIC_BC:
        case TETRAGONAL_BC:

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
            rmg::error("Lattice type not programmed");

    }                           /* end switch */


}                               /* end mg_restrict */




template <typename RmgType>
void mgrid::mg_prolong (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
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
        case -CUBIC_BC:
        case CUBIC_BC:
        case ORTHORHOMBIC_PRIMITIVE:
        case MONOCLINIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:
        case TETRAGONAL_BC:
        case TRICLINIC_PRIMITIVE:
        case No_Lattice:
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

        case HEXAGONAL2:
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
            rmg::error("Lattice type not programmed");

    } // end switch

}                               /* end mg_prolong */


template <typename RmgType>
void mgrid::mg_prolong_cubic (RmgType * __restrict__ full, RmgType * __restrict__ half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
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


    //rmg::pack_stop (half, rptr, dx2, dy2, dz2);
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
void mgrid::eval_residual (RmgType * __restrict__ mat, 
                           RmgType * __restrict__ f_mat, 
                           RmgType * __restrict__ work, 
                           int dimx, int dimy, int dimz,
                           double gridhx, double gridhy, double gridhz, 
                           RmgType * res, double *pot)
{
    int size, idx;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    FD.app_combined<RmgType, 2> (mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec.data(), false);
    rmg::pack_ptos(res, work, dimx, dimy, dimz);
    if(pot) {
        for (idx = 0; idx < size; idx++) res[idx] = f_mat[idx] + (RmgType)pot[idx]*mat[idx] - res[idx];
    }
    else {
        for (idx = 0; idx < size; idx++) res[idx] = f_mat[idx] - res[idx];
    }


}                               /* end eval_residual */



template <typename RmgType>
void mgrid::solv_pois (RmgType * __restrict__ vmat, RmgType * __restrict__ fmat, RmgType * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double k, double *pot)
{
    int size, idx;
    double scale;
    double diag;
    FiniteDiff FD(L);

    size = (dimx + 2) * (dimy + 2) * (dimz + 2);
    RmgType *work1 = &work[size];
//for (idx = 0; idx < size; idx++) work1[idx] = 0.0;

    diag = -FD.app_combined<RmgType, 2> (vmat, work1, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec.data(), false);
    rmg::pack_ptos(work, work1, dimx, dimy, dimz);
    double Zfac = zmax*gridhx*gridhx;
    scale = 1.0 / (diag + Zfac);
//printf("SCALE = %14.8f  %14.8f\n",step,scale);
    scale = step * scale;
//scale = 0.66*scale;
 
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

int mgrid::MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
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

    rmg::error("Boundary condition not programmed."); 
    return -1;

}

}

