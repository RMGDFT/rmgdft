/*
 *
 * Copyright (c) 1995,2011,2014 Emil Briggs
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

#include <cmath>
#include <complex>
#include <type_traits>
#include "Lattice.h"
#include "FiniteDiff.h"
#include "Gpufuncs.h"
#include "RmgTimer.h"
#include "rmg_error.h"

#define         PI          3.14159265358979323

int FiniteDiff::allocation_limit = 65536;
double FiniteDiff::cfac[12];
std::unordered_map<int, LaplacianCoeff *> FiniteDiff::FdCoeffs;

// Force instantiation of float, double and complex versions.
template double FiniteDiff::app_del2_np<double>(double *, double *, double, double, double);

template double FiniteDiff::app2_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template double FiniteDiff::app2_del2_offset<float>(float *, float *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<double>(double *, double *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

FiniteDiff::FiniteDiff(Lattice *lptr)
{
    L = lptr;
    this->xoff = NULL;
    this->yoff = NULL;
    this->zoff = NULL;
    this->x_type = PERIODIC;
    this->y_type = PERIODIC;
    this->z_type = PERIODIC;
    this->np_weight = NULL;
    this->alt_laplacian = false;
}

FiniteDiff::FiniteDiff(Lattice *lptr, bool alt_flag)
{
    L = lptr;
    this->xoff = NULL;
    this->yoff = NULL;
    this->zoff = NULL;
    this->np_weight = NULL;
    this->x_type = PERIODIC;
    this->y_type = PERIODIC;
    this->z_type = PERIODIC;
    this->alt_laplacian = alt_flag;
}

void FiniteDiff::set_allocation_limit(int lim)
{
    FiniteDiff::allocation_limit = lim;
}

// Constructor for non-periodic boundary conditions. Unlike the case with
// the standard constructor the non-periodic case is specific to a specific
// grid density and order of accuracy.
FiniteDiff::FiniteDiff(Lattice *lptr, BaseGrid *gptr, int xtype, int ytype, int ztype, int density, int norder)
{
    L = lptr;
    G = gptr;
    this->np_weight = NULL;
    this->xoff = NULL;
    this->yoff = NULL;
    this->zoff = NULL;
    this->np_density = density;
    this->x_type = PERIODIC;
    this->y_type = PERIODIC;
    this->z_type = PERIODIC;
    this->order = norder;

    // The stride var is the number of weights for each grid point. 
    // An eighth order second derivative for example requires 9 points in
    // the symmetric case and 10 points for the non-symmetric. We use 10 for
    // both in order to simplify later code and pad with a zero for the
    // symmetric case.
    this->stride = order + 2;
    if((xtype == PERIODIC) && (ytype == PERIODIC) && (ztype == PERIODIC)) return;

    this->np_weight = new double[2*order + 1]();

    this->xoff = new int[G->get_NX_GRID(this->np_density)];
    this->yoff = new int[G->get_NY_GRID(this->np_density)];
    this->zoff = new int[G->get_NZ_GRID(this->np_density)];

    int alloc = std::max(G->get_NX_GRID(this->np_density), G->get_NY_GRID(this->np_density));
    alloc = std::max(alloc, G->get_NZ_GRID(this->np_density));
    double *xr = new double[alloc];
    double *tw = new double[alloc];
    for(int i=0;i < alloc;i++) xr[i] = (double)i;

    int range = order / 2;

    // Default is periodic boundary condition
    for(int i =0 ;i < G->get_NX_GRID(this->np_density);i++) xoff[i] = -range;
    for(int i =0 ;i < G->get_NY_GRID(this->np_density);i++) yoff[i] = -range;
    for(int i =0 ;i < G->get_NZ_GRID(this->np_density);i++) zoff[i] = -range;
    for(int i=0;i < range;i++) xoff[i] = -i;
    for(int i=0;i < range;i++) yoff[i] = -i;
    for(int i=0;i < range;i++) zoff[i] = -i;

    //for(int i=0;i < order;i++) xoff[G->get_NX_GRID(this->np_density)-i-1] = xoff[i];
    //for(int i=0;i < order;i++) yoff[G->get_NY_GRID(this->np_density)-i-1] = yoff[i];
    //for(int i=0;i < order;i++) zoff[G->get_NZ_GRID(this->np_density)-i-1] = zoff[i];

    gen_weights(order + 1 , 2, (double)(order/2), &xr[0], tw);
    for(int j=0;j <= order;j++) this->np_weight[j] = tw[j];

#if 1
    //printf("\n");
    for(int j=0;j<=order;j++)
    {
        //printf("PPPPP    %d   %d   %d   %18.12f\n",j,order,this->xoff[j],this->np_weight[j]);
    }
#endif
    delete [] tw;
    delete [] xr;
}

//
void FiniteDiff::gen_weights(int n, int m, double xr, double *x, double *w)
{
    double c1 = 1.0, c2, c3, c4, c5;
    double c[20][20];

    c4 = x[0] - xr;

    for(int k=0;k <= m;k++)
    {
        for(int j=0;j <= n;j++)
        {
            c[j][k] = 0.0;
        }
    }
    c[0][0] = 1.0;

    for(int i=1;i <= n;i++)
    {
        int mn = std::min(i,m);
        c2 = 1.0;
        c5 = c4;

        c4 = x[i] - xr;

        for(int j=0;j <= i-1;j++)
        {
            c3 = x[i] - x[j];
            c2 = c2*c3;

            if (j == (i-1))
            {
                for(int k=mn;k >= 1;k--)
                {
                    c[i][k] = c1*(k*c[i-1][k-1] - c5*c[i-1][k]) / c2;
                }
                c[i][0] = -c1*c5*c[i-1][0] / c2;
            }
            for(int k=mn;k >= 1;k--)
            {
                c[j][k] = (c4*c[j][k]-k*c[j][k-1]) / c3;
            }
            c[j][0] = c4*c[j][0] / c3;
        }
        c1 = c2;
     }

     for(int i=0;i <= n;i++)w[i] = c[i][2];
}

void FiniteDiff::set_alt_laplacian_flag(bool flag)
{
    this->alt_laplacian = flag;
}

FiniteDiff::~FiniteDiff(void)
{
    if(this->xoff) delete [] this->xoff;
    if(this->yoff) delete [] this->yoff;
    if(this->zoff) delete [] this->zoff;
    if(this->np_weight) delete [] this->np_weight;
}

bool FiniteDiff::check_anisotropy(double hx, double hy, double hz, double limit)
{
    double anisotropy = hx / hy;
    if(fabs(anisotropy - 1.0) > limit) return false;
    anisotropy = hy / hz;
    if(fabs(anisotropy - 1.0) > limit) return false;
    anisotropy = hz / hx;
    if(fabs(anisotropy - 1.0) > limit) return false;
    return true;
}


// For the non-periodic case grid dims and spacing are set at construction time
template <typename RmgType>
double FiniteDiff::app_del2_np (RmgType *rptr, RmgType *b, double gridhx, double gridhy, double gridhz)
{
    int xdim = this->G->get_NX_GRID(this->np_density);
    int ydim = this->G->get_NY_GRID(this->np_density);
    int zdim = this->G->get_NZ_GRID(this->np_density);

    double ihx = 1.0 / (gridhx * this->L->get_xside());
    double ihy = 1.0 / (gridhy * this->L->get_yside());
    double ihz = 1.0 / (gridhz * this->L->get_zside());

    ihx = ihx*ihx;
    ihy = ihy*ihy;
    ihz = ihz*ihz;

    int incy = zdim;
    int incx = ydim*zdim;
    int range = order/ 2;

    for(int ix = 0;ix < xdim;ix++)
    {
        int xstart = 0;
        if(ix < range) xstart = range - ix;
        int xstop = order;
        if(ix > (xdim - range - 1)) xstop = range + (xdim - ix - 1);
//printf("XSTART  %d  %d  %d  %d\n",ix, xstart, xstop,xoff[ix]);
        for(int iy = 0;iy < ydim;iy++)
        {
            int ystart = 0;
            if(iy < range) ystart = range - iy;
            int ystop = order;
            if(iy > (ydim - range - 1)) ystop = range + (ydim - iy - 1);
            for(int iz = 0;iz < zdim;iz++)
            {
                int zstart = 0;
                if(iz < range) zstart = range - iz;
                int zstop = order;
                if(iz > (zdim - range - 1)) zstop = range + (zdim - iz - 1);
                int idx = ix*incx + iy*incy + iz;
                RmgType sumx = 0.0;
                RmgType sumy = 0.0;
                RmgType sumz = 0.0;
                for(int ii = xstart;ii <= xstop;ii++) sumx += this->np_weight[ii] * rptr[incx*(ix + ii + this->xoff[ix]) + iy*incy + iz];
                for(int ii = ystart;ii <= ystop;ii++) sumy += this->np_weight[ii] * rptr[ix*incx + incy*(iy + ii + this->yoff[iy]) + iz];
                for(int ii = zstart;ii <= zstop;ii++) sumz += this->np_weight[ii] * rptr[ix*incx + iy*incy + this->zoff[iz] + ii + iz];
                b[idx] = ihx*sumx + ihy*sumy + ihz*sumz;

            }
        }
    }
    double t0 = this->np_weight[this->order/2]*(ihx + ihy + ihz);
    return t0;
}


template <typename RmgType>
double FiniteDiff::app2_del2 (RmgType * __restrict__ a, RmgType * __restrict__  b, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz)
{

    int ix, iy, iz, ibrav;
    int incy, incx;
    RmgType cc=0.0, fcx, fcy, fcz, fc, fc1, fc2;
    int ixs, iys, ixms, ixps, iyms, iyps;
    RmgType ihx, ihy, ihz;
    RmgType ONE_t(1.0);
    RmgType TWO_t(2.0);
    RmgType FOUR_t(4.0);
    RmgType SIX_t(6.0);

    ihx = 1.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    ihy = 1.0 / (gridhy * gridhy * L->get_yside() * L->get_yside());
    ihz = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());

    ibrav = L->get_ibrav_type();


    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    switch (ibrav) {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:
        case TRICLINIC_PRIMITIVE:
        case MONOCLINIC_PRIMITIVE:
        case None:
            if (FiniteDiff::check_anisotropy(gridhx, gridhy, gridhz, 0.0000001))
            {

                cc = - TWO_t * ihx;
                cc = cc - TWO_t * ihy;
                cc = cc - TWO_t * ihz;
                fcx = ONE_t * ihx;

                for (ix = 1; ix <= dimx; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;

                    if (dimy % 2)
                    {

                        for (iy = 1; iy <= dimy; iy++)
                        {
                            iys = iy * incy;
                            iyms = (iy - 1) * incy;
                            iyps = (iy + 1) * incy;

                            for (iz = 1; iz <= dimz; iz++)
                            {

                                b[ixs + iys + iz] =
                                    cc * a[ixs + iys + iz] +
                                    fcx * (a[ixms + iys + iz] +
                                           a[ixps + iys + iz] +
                                           a[ixs + iyms + iz] +
                                           a[ixs + iyps + iz] +
                                           a[ixs + iys + (iz - 1)] + a[ixs + iys + (iz + 1)]);


                            }       /* end for */

                        }           /* end for */

                    }
                    else
                    {

                        for (iy = 1; iy <= dimy; iy += 2)
                        {
                            iys = iy * incy;
                            iyms = (iy - 1) * incy;
                            iyps = (iy + 1) * incy;

                            for (iz = 1; iz <= dimz; iz++)
                            {

                                b[ixs + iys + iz] =
                                    cc * a[ixs + iys + iz] +
                                    fcx * (a[ixs + iys + (iz - 1)] +
                                           a[ixs + iys + (iz + 1)] +
                                           a[ixms + iys + iz] +
                                           a[ixps + iys + iz] +
                                           a[ixs + iyms + iz] + a[ixs + iyps + iz]);

                                b[ixs + iyps + iz] =
                                    cc * a[ixs + iyps + iz] +
                                    fcx * (a[ixs + iyps + (iz - 1)] +
                                           a[ixs + iyps + (iz + 1)] +
                                           a[ixms + iyps + iz] +
                                           a[ixps + iyps + iz] +
                                           a[ixs + iys + iz] + a[ixs + iyps + incy + iz]);

                            }       /* end for */

                        }           /* end for */

                    }               /* end if */

                }                   /* end for */

            }
            else
            {

                cc = - TWO_t * ihx;
                cc = cc - TWO_t * ihy;
                cc = cc - TWO_t * ihz;
                fcx = ONE_t * ihx;
                fcy = ONE_t * ihy;
                fcz = ONE_t * ihz;

                for (ix = 1; ix <= dimx; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;

                    for (iy = 1; iy <= dimy; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (iz = 1; iz <= dimz; iz++)
                        {

                            b[ixs + iys + iz] =
                                cc * a[ixs + iys + iz] +
                                fcx * a[ixms + iys + iz] +
                                fcx * a[ixps + iys + iz] +
                                fcy * a[ixs + iyms + iz] +
                                fcy * a[ixs + iyps + iz] +
                                fcz * a[ixs + iys + (iz - 1)] + fcz * a[ixs + iys + (iz + 1)];


                        }           /* end for */

                    }               /* end for */

                }                   /* end for */

            }                       /* end if */

            break;

        case CUBIC_BC:

            cc = -TWO_t * ihx;
            fc = ONE_t * ihx / FOUR_t;

            for (ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (iz = 1; iz <= dimz; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc * a[ixms + iyms + iz - 1] +
                            fc * a[ixms + iys + iz] +
                            fc * a[ixs + iyms + iz] +
                            fc * a[ixs + iys + iz - 1] +
                            fc * a[ixs + iys + iz + 1] +
                            fc * a[ixs + iyps + iz] +
                            fc * a[ixps + iys + iz] + fc * a[ixps + iyps + iz + 1];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        case CUBIC_FC:

            cc = - SIX_t * ihx;
            fc = (ONE_t / TWO_t) * ihx; 


            for (ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (iz = 1; iz <= dimz; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc * a[ixms + iys + iz] +
                            fc * a[ixms + iys + iz + 1] +
                            fc * a[ixms + iyps + iz] +
                            fc * a[ixs + iyms + iz] +
                            fc * a[ixs + iyms + iz + 1] +
                            fc * a[ixs + iys + iz - 1] +
                            fc * a[ixs + iys + iz + 1] +
                            fc * a[ixs + iyps + iz - 1] +
                            fc * a[ixs + iyps + iz] +
                            fc * a[ixps + iyms + iz] +
                            fc * a[ixps + iys + iz - 1] + 
                            fc * a[ixps + iys + iz];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            cc = -FOUR_t * ihx;
            cc = cc - TWO_t * ihz;
            fc1 = 2.0 / (3.0 * gridhx * gridhx * L->get_xside() * L->get_xside());
            fc2 = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());

            for (ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (iz = 1; iz <= dimz; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc1 * a[ixps + iys + iz] +
                            fc1 * a[ixps + iyms + iz] +
                            fc1 * a[ixs + iyms + iz] +
                            fc1 * a[ixms + iys + iz] +
                            fc1 * a[ixms + iyps + iz] +
                            fc1 * a[ixs + iyps + iz] +
                            fc2 * a[ixs + iys + iz + 1] + fc2 * a[ixs + iys + iz - 1];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }                           /* end switch */


    /* Return the diagonal component of the operator */
    return (double)std::real(cc);

}                               /* end app2_del2 */


template <typename RmgType>
double FiniteDiff::app2_del2_offset (RmgType * a, RmgType * b, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz, int offset)
{

    int ibrav;
    int incy, incx;
    RmgType cc=0.0, fcx, fcy, fcz, fc, fc1, fc2;
    int ixs, iys, ixms, ixps, iyms, iyps;
    RmgType ihx, ihy, ihz;
    RmgType ONE_t(1.0);
    RmgType TWO_t(2.0);
    RmgType FOUR_t(4.0);
    RmgType SIX_t(6.0);

    ihx = 1.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    ihy = 1.0 / (gridhy * gridhy * L->get_yside() * L->get_yside());
    ihz = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());

    ibrav = L->get_ibrav_type();


    incy = dimz;
    incx = dimz * dimy;

    switch (ibrav) {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:
        case TRICLINIC_PRIMITIVE:

            if (FiniteDiff::check_anisotropy(gridhx, gridhy, gridhz, 0.0000001))
            {

                cc = - TWO_t * ihx;
                cc = cc - TWO_t * ihy;
                cc = cc - TWO_t * ihz;
                fcx = ONE_t * ihx;

                for (int ix = offset; ix < dimx-offset; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;

                    for (int iy = offset; iy < dimy-offset; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (int iz = offset; iz < dimz-offset; iz++)
                        {

                            b[ixs+iys+iz] =
                                cc * a[ixs + iys + iz] +
                                fcx * (a[ixms + iys + iz] +
                                       a[ixps + iys + iz] +
                                       a[ixs + iyms + iz] +
                                       a[ixs + iyps + iz] +
                                       a[ixs + iys + (iz - 1)] + 
                                       a[ixs + iys + (iz + 1)]);

                        }       /* end for */

                    }           /* end for */

                }                   /* end for */

            }
            else
            {

                cc = - TWO_t * ihx;
                cc = cc - TWO_t * ihy;
                cc = cc - TWO_t * ihz;
                fcx = ONE_t * ihx;
                fcy = ONE_t * ihy;
                fcz = ONE_t * ihz;

                for (int ix = offset; ix < dimx-offset; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;

                    for (int iy = offset; iy < dimy-offset; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (int iz = offset; iz < dimz-offset; iz++)
                        {

                            b[ixs + iys + iz] =
                                cc * a[ixs + iys + iz] +
                                fcx * a[ixms + iys + iz] +
                                fcx * a[ixps + iys + iz] +
                                fcy * a[ixs + iyms + iz] +
                                fcy * a[ixs + iyps + iz] +
                                fcz * a[ixs + iys + (iz - 1)] + 
                                fcz * a[ixs + iys + (iz + 1)];


                        }           /* end for */

                    }               /* end for */

                }                   /* end for */

            }                       /* end if */

            break;

        case CUBIC_BC:

            cc = -TWO_t * ihx;
            fc = ONE_t * ihx / FOUR_t;

            for (int ix = offset; ix < dimx-offset; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (int iy = offset; iy < dimy-offset; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = offset; iz < dimz-offset; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc * a[ixms + iyms + iz - 1] +
                            fc * a[ixms + iys + iz] +
                            fc * a[ixs + iyms + iz] +
                            fc * a[ixs + iys + iz - 1] +
                            fc * a[ixs + iys + iz + 1] +
                            fc * a[ixs + iyps + iz] +
                            fc * a[ixps + iys + iz] + fc * a[ixps + iyps + iz + 1];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        case CUBIC_FC:

            cc = - SIX_t * ihx;
            fc = (ONE_t / TWO_t) * ihx; 

            for (int ix = offset; ix < dimx-offset; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (int iy = offset; iy < dimy-offset; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = offset; iz < dimz-offset; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc * a[ixms + iys + iz] +
                            fc * a[ixms + iys + iz + 1] +
                            fc * a[ixms + iyps + iz] +
                            fc * a[ixs + iyms + iz] +
                            fc * a[ixs + iyms + iz + 1] +
                            fc * a[ixs + iys + iz - 1] +
                            fc * a[ixs + iys + iz + 1] +
                            fc * a[ixs + iyps + iz - 1] +
                            fc * a[ixs + iyps + iz] +
                            fc * a[ixps + iyms + iz] +
                            fc * a[ixps + iys + iz - 1] + 
                            fc * a[ixps + iys + iz];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            cc = -FOUR_t * ihx;
            cc = cc - TWO_t * ihz;
            fc1 = 2.0 / (3.0 * gridhx * gridhx * L->get_xside() * L->get_xside());
            fc2 = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());

            for (int ix = offset; ix < dimx-offset; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (int iy = offset; iy < dimy-offset; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = offset; iz < dimz-offset; iz++)
                    {

                        b[ixs + iys + iz] =
                            cc * a[ixs + iys + iz] +
                            fc1 * a[ixps + iys + iz] +
                            fc1 * a[ixps + iyms + iz] +
                            fc1 * a[ixs + iyms + iz] +
                            fc1 * a[ixms + iys + iz] +
                            fc1 * a[ixms + iyps + iz] +
                            fc1 * a[ixs + iyps + iz] +
                            fc2 * a[ixs + iys + iz + 1] + fc2 * a[ixs + iys + iz - 1];


                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }                           /* end switch */


    /* Return the diagonal component of the operator */
    return (double)std::real(cc);

}                               /* end app2_del2_offset */

