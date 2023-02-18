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
double FiniteDiff::cfac[13];
std::unordered_map<int, LaplacianCoeff *> FiniteDiff::FdCoeffs;


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

// Generates a key for the FdCoeffs map
// a0h is the grid spacing for the first axis
int FiniteDiff::LCkey(double a0h)
{
    return (int)std::round(a0h * 100000);
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


#include "rmg_complex.h"


// Force instantiation of float, double and complex versions.
template double FiniteDiff::app2_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template double FiniteDiff::app8_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template void FiniteDiff::app_gradient_eighth<float> (float *, float *, float *, float *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_eighth<double> (double *, double *, double *, double *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_eighth<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int, double , double , double );
template void FiniteDiff::app_gradient_eighth<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int, double , double , double );

template void FiniteDiff::fd_gradient_general<float, 6> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<double, 6> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<double>, 6> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<float>, 6> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general<float, 8> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<double, 8> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<double>, 8> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<float>, 8> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general<float, 10> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<double, 10> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<double>, 10> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<float>, 10> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general<float, 12> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<double, 12> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<double>, 12> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general<std::complex<float>, 12> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template double FiniteDiff::app_combined<float,2>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,2>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 2>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 2>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<float,4>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,4>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 4>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 4>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<float,6>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,6>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 6>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 6>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<float,8>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,8>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 8>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 8>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<float,10>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,10>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 10>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 10>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

template double FiniteDiff::app_combined<float,12>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<double,12>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <float>, 12>(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app_combined<std::complex <double>, 12>(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

template void FiniteDiff::fd_gradient_coeffs<float>(int , double, int , float *, float *, float *);
template void FiniteDiff::fd_gradient_coeffs<double>(int , double, int , double *, double *, double *);
template void FiniteDiff::fd_gradient_coeffs<std::complex<double>>(int , double, int , std::complex<double> *, std::complex<double> *, std::complex<double> *);

template void FiniteDiff::fd_combined_coeffs<float>(int, double, int, float *, float *, double *);
template void FiniteDiff::fd_combined_coeffs<double>(int, double, int, double *, double *, double *);
template void FiniteDiff::fd_combined_coeffs<std::complex<float>>(int, double, int, std::complex<float> *, std::complex<float> *, double *);
template void FiniteDiff::fd_combined_coeffs<std::complex<double>>(int, double, int, std::complex<double> *, std::complex<double> *, double *);


template <typename RmgType>
double FiniteDiff::app2_del2 (RmgType * __restrict__ a, RmgType * __restrict__  b, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz)
{
    double kvec[3] = {0.0, 0.0, 0.0};
    return FiniteDiff::app_combined<RmgType, 2>(a, b, dimx, dimy, dimz,
        gridhx, gridhy, gridhz, kvec, false);
}                               /* end app2_del2 */


template <typename RmgType>
double FiniteDiff::app8_del2(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
               double gridhx, double gridhy, double gridhz)
{
    double kvec[3] = {0.0,0.0,0.0};
    return FiniteDiff::app_combined<RmgType, 8>(a, b, dimx, dimy, dimz,
        gridhx, gridhy, gridhz, kvec, false);

}  /* end app8_del2 */


template <typename RmgType>
void FiniteDiff::app_gradient_eighth (RmgType * __restrict__ rptr, RmgType * __restrict__ wxr, RmgType * __restrict__ wyr, RmgType * __restrict__ wzr, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{
    FiniteDiff::fd_gradient_general<RmgType, 8> (rptr, wxr, wyr, wzr, gridhx, dimx, dimy, dimz);
}


template <typename RmgType, int order>
double FiniteDiff::app_combined(RmgType * __restrict__ a, RmgType * __restrict__ b, 
		int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
		double *kvec, bool use_gpu)
{
    int ibrav = L->get_ibrav_type();
    RmgType cm[12], cp[12];
    RmgType cpx[12], cmx[12], cpy[12], cmy[12], cpz[12], cmz[12];
    int ixs = (dimy + order) * (dimz + order);
    int iys = (dimz + order);


    // NULL b means we just want the diagonal component.
    double th2 = fd_coeff0(order, gridhx);
    if(b == NULL) return (double)std::real(th2);

#if 0
#if HIP_ENABLED || CUDA_ENABLED
    // Broken for now. Need to set up c
    if(use_gpu && (ibrav == CUBIC_PRIMITIVE || ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE))
    {
        /* Return the diagonal component of the operator */
        app8_del2_gpu(a, b, dimx, dimy, dimz, c);
        return (double)std::real(th2);
    }
#endif
#endif

    // Get coeffs for x,y,z axes which are used by all lattice types
    fd_combined_coeffs(order, gridhx, 0, cmx, cpx, kvec);
    fd_combined_coeffs(order, gridhx, 1, cmy, cpy, kvec);
    fd_combined_coeffs(order, gridhx, 2, cmz, cpz, kvec);

    for (int ix = order/2; ix < dimx + order/2; ix++)
    {
        for (int iy = order/2; iy < dimy + order/2; iy++)
        {
            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
            // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
            // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                B[iz] = th2 * A[iz];
                B[iz] += cpz[0] * A[iz + 1] + cmz[0] * A[iz - 1];
                if constexpr(order >= 4)
                    B[iz] += cpz[1] * A[iz + 2] + cmz[1] * A[iz - 2];
                if constexpr(order >= 6)
                    B[iz] += cpz[2] * A[iz + 3] + cmz[2] * A[iz - 3];
                if constexpr(order >= 8)
                    B[iz] += cpz[3] * A[iz + 4] + cmz[3] * A[iz - 4];
                if constexpr(order >= 10)
                    B[iz] += cpz[4] * A[iz + 5] + cmz[4] * A[iz - 5];
                if constexpr(order >= 12)
                    B[iz] += cpz[5] * A[iz + 6] + cmz[5] * A[iz - 6];
            }
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                B[iz] += cpy[0] * A[iz + iys] + cmy[0] * A[iz - iys];
                if constexpr(order >= 4)
                    B[iz] += cpy[1] * A[iz + 2*iys] + cmy[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    B[iz] += cpy[2] * A[iz + 3*iys] + cmy[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    B[iz] += cpy[3] * A[iz + 4*iys] + cmy[3] * A[iz - 4*iys];
                if constexpr(order >= 10)
                    B[iz] += cpy[4] * A[iz + 5*iys] + cmy[4] * A[iz - 5*iys];
                if constexpr(order >= 12)
                    B[iz] += cpy[5] * A[iz + 6*iys] + cmy[5] * A[iz - 6*iys];
            }
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                B[iz] += cpx[0] * A[iz + ixs] + cmx[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    B[iz] += cpx[1] * A[iz + 2*ixs] + cmx[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    B[iz] += cpx[2] * A[iz + 3*ixs] + cmx[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    B[iz] += cpx[3] * A[iz + 4*ixs] + cmx[3] * A[iz - 4*ixs];
                if constexpr(order >= 10)
                    B[iz] += cpx[4] * A[iz + 5*ixs] + cmx[4] * A[iz - 5*ixs];
                if constexpr(order >= 12)
                    B[iz] += cpx[5] * A[iz + 6*ixs] + cmx[5] * A[iz - 6*ixs];
            }                   /* end for */
        }
    }

    /* Quick return for orthogonal axis cases */
    if(ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == CUBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
        return (double)std::real(th2);

    // Add additional axes as required
    if(LC->include_axis[3])
    {
        fd_combined_coeffs(order, gridhx, 3, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + ixs + iys] + cm[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*ixs + 2*iys] + cm[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*ixs + 3*iys] + cm[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*ixs + 4*iys] + cm[3] * A[iz - 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*ixs + 5*iys] + cm[4] * A[iz - 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*ixs + 6*iys] + cm[5] * A[iz - 6*ixs - 6*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[4])
    {
        fd_combined_coeffs(order, gridhx, 4, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + ixs + 1] + cm[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*ixs + 2] + cm[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*ixs + 3] + cm[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*ixs + 4] + cm[3] * A[iz - 4*ixs - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*ixs + 5] + cm[4] * A[iz - 5*ixs - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*ixs + 6] + cm[5] * A[iz - 6*ixs - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[5])
    {
        fd_combined_coeffs(order, gridhx, 5, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + iys + 1] + cm[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*iys + 2] + cm[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*iys + 3] + cm[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*iys + 4] + cm[3] * A[iz - 4*iys - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*iys + 5] + cm[4] * A[iz - 5*iys - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*iys + 6] + cm[5] * A[iz - 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[6])
    {
        fd_combined_coeffs(order, gridhx, 6, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz - ixs + iys] + cm[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz - 2*ixs + 2*iys] + cm[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz - 3*ixs + 3*iys] + cm[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz - 4*ixs + 4*iys] + cm[3] * A[iz + 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz - 5*ixs + 5*iys] + cm[4] * A[iz + 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz - 6*ixs + 6*iys] + cm[5] * A[iz + 6*ixs - 6*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[7])
    {
        fd_combined_coeffs(order, gridhx, 7, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz - ixs + 1] + cm[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz - 2*ixs + 2] + cm[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz - 3*ixs + 3] + cm[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz - 4*ixs + 4] + cm[3] * A[iz + 4*ixs - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz - 5*ixs + 5] + cm[4] * A[iz + 5*ixs - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz - 6*ixs + 6] + cm[5] * A[iz + 6*ixs - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[8])
    {
        fd_combined_coeffs(order, gridhx, 8, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz - iys + 1] + cm[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz - 2*iys + 2] + cm[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz - 3*iys + 3] + cm[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz - 4*iys + 4] + cm[3] * A[iz + 4*iys - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz - 5*iys + 5] + cm[4] * A[iz + 5*iys - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz - 6*iys + 6] + cm[5] * A[iz + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[9])
    {
        fd_combined_coeffs(order, gridhx, 9, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + 1*ixs + 1*iys + 1] + cm[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*ixs + 2*iys + 2] + cm[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*ixs + 3*iys + 3] + cm[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*ixs + 4*iys + 4] + cm[3] * A[iz - 4*ixs - 4*iys - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*ixs + 5*iys + 5] + cm[4] * A[iz - 5*ixs - 5*iys - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*ixs + 6*iys + 6] + cm[5] * A[iz - 6*ixs - 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[10])
    {
        fd_combined_coeffs(order, gridhx, 10, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz - 1*ixs - 1*iys + 1] + cm[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz - 2*ixs - 2*iys + 2] + cm[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz - 3*ixs - 3*iys + 3] + cm[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz - 4*ixs - 4*iys + 4] + cm[3] * A[iz + 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz - 5*ixs - 5*iys + 5] + cm[4] * A[iz + 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz - 6*ixs - 6*iys + 6] + cm[5] * A[iz + 6*ixs + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[11])
    {
        fd_combined_coeffs(order, gridhx, 11, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + 1*ixs - 1*iys + 1] + cm[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*ixs - 2*iys + 2] + cm[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*ixs - 3*iys + 3] + cm[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*ixs - 4*iys + 4] + cm[3] * A[iz - 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*ixs - 5*iys + 5] + cm[4] * A[iz - 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*ixs - 6*iys + 6] + cm[5] * A[iz - 6*ixs + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[12])
    {
        fd_combined_coeffs(order, gridhx, 12, cm, cp, kvec);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    B[iz] += cp[0] * A[iz + 1*ixs - 1*iys - 1] + cm[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        B[iz] += cp[1] * A[iz + 2*ixs - 2*iys - 2] + cm[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        B[iz] += cp[2] * A[iz + 3*ixs - 3*iys - 3] + cm[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        B[iz] += cp[3] * A[iz + 4*ixs - 4*iys - 4] + cm[3] * A[iz - 4*ixs + 4*iys + 4];
                    if constexpr(order >= 10)
                        B[iz] += cp[4] * A[iz + 5*ixs - 5*iys - 5] + cm[4] * A[iz - 5*ixs + 5*iys + 5];
                    if constexpr(order >= 12)
                        B[iz] += cp[5] * A[iz + 6*ixs - 6*iys - 6] + cm[5] * A[iz - 6*ixs + 6*iys + 6];
                }                   /* end for */
            }
        }
    }

    /* Return the diagonal component of the operator */
    return (double)std::real(th2);

} /* end app_combined */



template <typename RmgType, int order>
void FiniteDiff::fd_gradient_general (RmgType * __restrict__ a, 
                                RmgType * __restrict__ gx, 
                                RmgType * __restrict__ gy, 
                                RmgType * __restrict__ gz,
                                double gridhx,
                                int dimx, int dimy, int dimz)
{
    int ibrav = L->get_ibrav_type();
    RmgType cxx[12], cxy[12], cxz[12];
    RmgType cyx[12], cyy[12], cyz[12];
    RmgType czx[12], czy[12], czz[12];
    RmgType cx[12], cy[12], cz[12];
    int ixs = (dimy + order) * (dimz + order);
    int iys = (dimz + order);
    LaplacianCoeff *LC = FiniteDiff::FdCoeffs[FiniteDiff::LCkey(gridhx) + order];

    // Get coeffs for x,y,z axes which are used by all lattice types
    fd_gradient_coeffs(order, gridhx, 0, cxx, cxy, cxz);
    fd_gradient_coeffs(order, gridhx, 1, cyx, cyy, cyz);
    fd_gradient_coeffs(order, gridhx, 2, czx, czy, czz);

    for (int ix = order/2; ix < dimx + order/2; ix++)
    {
        for (int iy = order/2; iy < dimy + order/2; iy++)
        {
            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
            RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
            RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
            // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
            // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                bgx[iz] = -czx[0] * A[iz + 1] + czx[0] * A[iz - 1];
                if constexpr(order >= 4)
                    bgx[iz] += -czx[1] * A[iz + 2] + czx[1] * A[iz - 2];
                if constexpr(order >= 6)
                    bgx[iz] += -czx[2] * A[iz + 3] + czx[2] * A[iz - 3];
                if constexpr(order >= 8)
                    bgx[iz] += -czx[3] * A[iz + 4] + czx[3] * A[iz - 4];
                if constexpr(order >= 10)
                    bgx[iz] += -czx[4] * A[iz + 5] + czx[4] * A[iz - 5];
                if constexpr(order >= 12)
                    bgx[iz] += -czx[5] * A[iz + 6] + czx[5] * A[iz - 6];

                bgy[iz] = -czy[0] * A[iz + 1] + czy[0] * A[iz - 1];
                if constexpr(order >= 4)
                    bgy[iz] += -czy[1] * A[iz + 2] + czy[1] * A[iz - 2];
                if constexpr(order >= 6)
                    bgy[iz] += -czy[2] * A[iz + 3] + czy[2] * A[iz - 3];
                if constexpr(order >= 8)
                    bgy[iz] += -czy[3] * A[iz + 4] + czy[3] * A[iz - 4];
                if constexpr(order >= 10)
                    bgy[iz] += -czy[4] * A[iz + 5] + czy[4] * A[iz - 5];
                if constexpr(order >= 12)
                    bgy[iz] += -czy[5] * A[iz + 6] + czy[5] * A[iz - 6];

                bgz[iz] = -czz[0] * A[iz + 1] + czz[0] * A[iz - 1];
                if constexpr(order >= 4)
                    bgz[iz] += -czz[1] * A[iz + 2] + czz[1] * A[iz - 2];
                if constexpr(order >= 6)
                    bgz[iz] += -czz[2] * A[iz + 3] + czz[2] * A[iz - 3];
                if constexpr(order >= 8)
                    bgz[iz] += -czz[3] * A[iz + 4] + czz[3] * A[iz - 4];
                if constexpr(order >= 10)
                    bgz[iz] += -czz[4] * A[iz + 5] + czz[4] * A[iz - 5];
                if constexpr(order >= 12)
                    bgz[iz] += -czz[5] * A[iz + 6] + czz[5] * A[iz - 6];

            }
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                bgx[iz] += -cyx[0] * A[iz + iys] + cyx[0] * A[iz - iys];
                if constexpr(order >= 4)
                    bgx[iz] += -cyx[1] * A[iz + 2*iys] + cyx[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    bgx[iz] += -cyx[2] * A[iz + 3*iys] + cyx[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    bgx[iz] += -cyx[3] * A[iz + 4*iys] + cyx[3] * A[iz - 4*iys];
                if constexpr(order >= 10)
                    bgx[iz] += -cyx[4] * A[iz + 5*iys] + cyx[4] * A[iz - 5*iys];
                if constexpr(order >= 12)
                    bgx[iz] += -cyx[5] * A[iz + 6*iys] + cyx[5] * A[iz - 6*iys];

                bgy[iz] += -cyy[0] * A[iz + iys] + cyy[0] * A[iz - iys];
                if constexpr(order >= 4)
                    bgy[iz] += -cyy[1] * A[iz + 2*iys] + cyy[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    bgy[iz] += -cyy[2] * A[iz + 3*iys] + cyy[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    bgy[iz] += -cyy[3] * A[iz + 4*iys] + cyy[3] * A[iz - 4*iys];
                if constexpr(order >= 10)
                    bgy[iz] += -cyy[4] * A[iz + 5*iys] + cyy[4] * A[iz - 5*iys];
                if constexpr(order >= 12)
                    bgy[iz] += -cyy[5] * A[iz + 6*iys] + cyy[5] * A[iz - 6*iys];

                bgz[iz] += -cyz[0] * A[iz + iys] + cyz[0] * A[iz - iys];
                if constexpr(order >= 4)
                    bgz[iz] += -cyz[1] * A[iz + 2*iys] + cyz[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    bgz[iz] += -cyz[2] * A[iz + 3*iys] + cyz[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    bgz[iz] += -cyz[3] * A[iz + 4*iys] + cyz[3] * A[iz - 4*iys];
                if constexpr(order >= 10)
                    bgz[iz] += -cyz[4] * A[iz + 5*iys] + cyz[4] * A[iz - 5*iys];
                if constexpr(order >= 12)
                    bgz[iz] += -cyz[5] * A[iz + 6*iys] + cyz[5] * A[iz - 6*iys];
            }
            for (int iz = order/2; iz < dimz + order/2; iz++)
            {
                bgx[iz] += -cxx[0] * A[iz + ixs] + cxx[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    bgx[iz] += -cxx[1] * A[iz + 2*ixs] + cxx[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    bgx[iz] += -cxx[2] * A[iz + 3*ixs] + cxx[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    bgx[iz] += -cxx[3] * A[iz + 4*ixs] + cxx[3] * A[iz - 4*ixs];
                if constexpr(order >= 10)
                    bgx[iz] += -cxx[4] * A[iz + 5*ixs] + cxx[4] * A[iz - 5*ixs];
                if constexpr(order >= 12)
                    bgx[iz] += -cxx[5] * A[iz + 6*ixs] + cxx[5] * A[iz - 6*ixs];

                bgy[iz] += -cxy[0] * A[iz + ixs] + cxy[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    bgy[iz] +=-cxy[1] * A[iz + 2*ixs] + cxy[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    bgy[iz] +=-cxy[2] * A[iz + 3*ixs] + cxy[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    bgy[iz] +=-cxy[3] * A[iz + 4*ixs] + cxy[3] * A[iz - 4*ixs];
                if constexpr(order >= 10)
                    bgy[iz] +=-cxy[4] * A[iz + 5*ixs] + cxy[4] * A[iz - 5*ixs];
                if constexpr(order >= 12)
                    bgy[iz] +=-cxy[5] * A[iz + 6*ixs] + cxy[5] * A[iz - 6*ixs];

                bgz[iz] += -cxz[0] * A[iz + ixs] + cxz[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    bgz[iz] += -cxz[1] * A[iz + 2*ixs] + cxz[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    bgz[iz] += -cxz[2] * A[iz + 3*ixs] + cxz[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    bgz[iz] += -cxz[3] * A[iz + 4*ixs] + cxz[3] * A[iz - 4*ixs];
                if constexpr(order >= 10)
                    bgz[iz] += -cxz[4] * A[iz + 5*ixs] + cxz[4] * A[iz - 5*ixs];
                if constexpr(order >= 12)
                    bgz[iz] += -cxz[5] * A[iz + 6*ixs] + cxz[5] * A[iz - 6*ixs];

            }                   /* end for */
        }
    }

    /* Quick return for orthogonal axis cases */
    if(ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == CUBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
        return;

    // Add additional axes as required
    if(LC->include_axis[3])
    {
        fd_gradient_coeffs(order, gridhx, 3, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + ixs + iys] + cz[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*ixs + 2*iys] + cz[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*ixs + 3*iys] + cz[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*ixs + 4*iys] + cz[3] * A[iz - 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*ixs + 5*iys] + cz[4] * A[iz - 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*ixs + 6*iys] + cz[5] * A[iz - 6*ixs - 6*iys];

                    bgy[iz] += -cy[0] * A[iz + ixs + iys] + cy[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2*iys] + cy[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3*iys] + cy[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4*iys] + cy[3] * A[iz - 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*ixs + 5*iys] + cy[4] * A[iz - 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*ixs + 6*iys] + cy[5] * A[iz - 6*ixs - 6*iys];

                    bgx[iz] += -cx[0] * A[iz + ixs + iys] + cx[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2*iys] + cx[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3*iys] + cx[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4*iys] + cx[3] * A[iz - 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*ixs + 5*iys] + cx[4] * A[iz - 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*ixs + 6*iys] + cx[5] * A[iz - 6*ixs - 6*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[4])
    {
        fd_gradient_coeffs(order, gridhx, 4, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + ixs + 1] + cz[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*ixs + 2] + cz[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*ixs + 3] + cz[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*ixs + 4] + cz[3] * A[iz - 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*ixs + 5] + cz[4] * A[iz - 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*ixs + 6] + cz[5] * A[iz - 6*ixs - 6];

                    bgy[iz] += -cy[0] * A[iz + ixs + 1] + cy[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2] + cy[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3] + cy[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4] + cy[3] * A[iz - 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*ixs + 5] + cy[4] * A[iz - 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*ixs + 6] + cy[5] * A[iz - 6*ixs - 6];

                    bgx[iz] += -cx[0] * A[iz + ixs + 1] + cx[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2] + cx[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3] + cx[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4] + cx[3] * A[iz - 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*ixs + 5] + cx[4] * A[iz - 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*ixs + 6] + cx[5] * A[iz - 6*ixs - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[5])
    {
        fd_gradient_coeffs(order, gridhx, 5, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + iys + 1] + cz[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*iys + 2] + cz[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*iys + 3] + cz[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*iys + 4] + cz[3] * A[iz - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*iys + 5] + cz[4] * A[iz - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*iys + 6] + cz[5] * A[iz - 6*iys - 6];

                    bgy[iz] += -cy[0] * A[iz + iys + 1] + cy[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*iys + 2] + cy[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*iys + 3] + cy[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*iys + 4] + cy[3] * A[iz - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*iys + 5] + cy[4] * A[iz - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*iys + 6] + cy[5] * A[iz - 6*iys - 6];

                    bgx[iz] += -cx[0] * A[iz + iys + 1] + cx[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*iys + 2] + cx[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*iys + 3] + cx[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*iys + 4] + cx[3] * A[iz - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*iys + 5] + cx[4] * A[iz - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*iys + 6] + cx[5] * A[iz - 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[6])
    {
        fd_gradient_coeffs(order, gridhx, 6, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz - ixs + iys] + cz[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz - 2*ixs + 2*iys] + cz[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz - 3*ixs + 3*iys] + cz[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz - 4*ixs + 4*iys] + cz[3] * A[iz + 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz - 5*ixs + 5*iys] + cz[4] * A[iz + 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*ixs + 6*iys] + cz[5] * A[iz + 6*ixs - 6*iys];

                    bgy[iz] += -cy[0] * A[iz - ixs + iys] + cy[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs + 2*iys] + cy[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs + 3*iys] + cy[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs + 4*iys] + cy[3] * A[iz + 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz - 5*ixs + 5*iys] + cy[4] * A[iz + 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*ixs + 6*iys] + cz[5] * A[iz + 6*ixs - 6*iys];

                    bgx[iz] += -cx[0] * A[iz - ixs + iys] + cx[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs + 2*iys] + cx[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs + 3*iys] + cx[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs + 4*iys] + cx[3] * A[iz + 4*ixs - 4*iys];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz - 5*ixs + 5*iys] + cx[4] * A[iz + 5*ixs - 5*iys];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*ixs + 6*iys] + cz[5] * A[iz + 6*ixs - 6*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[7])
    {
        fd_gradient_coeffs(order, gridhx, 7, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz - ixs + 1] + cz[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz - 2*ixs + 2] + cz[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz - 3*ixs + 3] + cz[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz - 4*ixs + 4] + cz[3] * A[iz + 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz - 5*ixs + 5] + cz[4] * A[iz + 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*ixs + 6] + cz[5] * A[iz + 6*ixs - 6];

                    bgy[iz] += -cy[0] * A[iz - ixs + 1] + cy[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs + 2] + cy[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs + 3] + cy[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs + 4] + cy[3] * A[iz + 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz - 5*ixs + 5] + cy[4] * A[iz + 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz - 6*ixs + 6] + cy[5] * A[iz + 6*ixs - 6];

                    bgx[iz] += -cx[0] * A[iz - ixs + 1] + cx[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs + 2] + cx[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs + 3] + cx[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs + 4] + cx[3] * A[iz + 4*ixs - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz - 5*ixs + 5] + cx[4] * A[iz + 5*ixs - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz - 6*ixs + 6] + cx[5] * A[iz + 6*ixs - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[8])
    {
        fd_gradient_coeffs(order, gridhx, 8, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz - iys + 1] + cz[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz - 2*iys + 2] + cz[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz - 3*iys + 3] + cz[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz - 4*iys + 4] + cz[3] * A[iz + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz - 5*iys + 5] + cz[4] * A[iz + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*iys + 6] + cz[5] * A[iz + 6*iys - 6];

                    bgy[iz] += -cy[0] * A[iz - iys + 1] + cy[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*iys + 2] + cy[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*iys + 3] + cy[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*iys + 4] + cy[3] * A[iz + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz - 5*iys + 5] + cy[4] * A[iz + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz - 6*iys + 6] + cy[5] * A[iz + 6*iys - 6];

                    bgx[iz] += -cx[0] * A[iz - iys + 1] + cx[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*iys + 2] + cx[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*iys + 3] + cx[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*iys + 4] + cx[3] * A[iz + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz - 5*iys + 5] + cx[4] * A[iz + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz - 6*iys + 6] + cx[5] * A[iz + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[9])
    {
        fd_gradient_coeffs(order, gridhx, 9, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + 1*ixs + 1*iys + 1] + cz[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*ixs + 2*iys + 2] + cz[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*ixs + 3*iys + 3] + cz[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*ixs + 4*iys + 4] + cz[3] * A[iz - 4*ixs - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*ixs + 5*iys + 5] + cz[4] * A[iz - 5*ixs - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*ixs + 6*iys + 6] + cz[5] * A[iz - 6*ixs - 6*iys - 6];

                    bgy[iz] += -cy[0] * A[iz + 1*ixs + 1*iys + 1] + cy[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2*iys + 2] + cy[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3*iys + 3] + cy[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4*iys + 4] + cy[3] * A[iz - 4*ixs - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*ixs + 5*iys + 5] + cy[4] * A[iz - 5*ixs - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*ixs + 6*iys + 6] + cy[5] * A[iz - 6*ixs - 6*iys - 6];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs + 1*iys + 1] + cx[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2*iys + 2] + cx[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3*iys + 3] + cx[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4*iys + 4] + cx[3] * A[iz - 4*ixs - 4*iys - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*ixs + 5*iys + 5] + cx[4] * A[iz - 5*ixs - 5*iys - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*ixs + 6*iys + 6] + cx[5] * A[iz - 6*ixs - 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[10])
    {
        fd_gradient_coeffs(order, gridhx, 10, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz - 1*ixs - 1*iys + 1] + cz[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz - 2*ixs - 2*iys + 2] + cz[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz - 3*ixs - 3*iys + 3] + cz[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz - 4*ixs - 4*iys + 4] + cz[3] * A[iz + 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz - 5*ixs - 5*iys + 5] + cz[4] * A[iz + 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz - 6*ixs - 6*iys + 6] + cz[5] * A[iz + 6*ixs + 6*iys - 6];

                    bgy[iz] += -cy[0] * A[iz - 1*ixs - 1*iys + 1] + cy[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs - 2*iys + 2] + cy[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs - 3*iys + 3] + cy[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs - 4*iys + 4] + cy[3] * A[iz + 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz - 5*ixs - 5*iys + 5] + cy[4] * A[iz + 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz - 6*ixs - 6*iys + 6] + cy[5] * A[iz + 6*ixs + 6*iys - 6];

                    bgx[iz] += -cx[0] * A[iz - 1*ixs - 1*iys + 1] + cx[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs - 2*iys + 2] + cx[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs - 3*iys + 3] + cx[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs - 4*iys + 4] + cx[3] * A[iz + 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz - 5*ixs - 5*iys + 5] + cx[4] * A[iz + 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz - 6*ixs - 6*iys + 6] + cx[5] * A[iz + 6*ixs + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[11])
    {
        fd_gradient_coeffs(order, gridhx, 11, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + 1*ixs - 1*iys + 1] + cz[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*ixs - 2*iys + 2] + cz[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*ixs - 3*iys + 3] + cz[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*ixs - 4*iys + 4] + cz[3] * A[iz - 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*ixs - 5*iys + 5] + cz[4] * A[iz - 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*ixs - 6*iys + 6] + cz[5] * A[iz - 6*ixs + 6*iys - 6];

                    bgy[iz] += -cy[0] * A[iz + 1*ixs - 1*iys + 1] + cy[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs - 2*iys + 2] + cy[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs - 3*iys + 3] + cy[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs - 4*iys + 4] + cy[3] * A[iz - 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*ixs - 5*iys + 5] + cy[4] * A[iz - 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*ixs - 6*iys + 6] + cy[5] * A[iz - 6*ixs + 6*iys - 6];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs - 1*iys + 1] + cx[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs - 2*iys + 2] + cx[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs - 3*iys + 3] + cx[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs - 4*iys + 4] + cx[3] * A[iz - 4*ixs + 4*iys - 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*ixs - 5*iys + 5] + cx[4] * A[iz - 5*ixs + 5*iys - 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*ixs - 6*iys + 6] + cx[5] * A[iz - 6*ixs + 6*iys - 6];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[12])
    {
        fd_gradient_coeffs(order, gridhx, 12, cx, cy, cz);
        for (int ix = order/2; ix < dimx + order/2; ix++)
        {
            for (int iy = order/2; iy < dimy + order/2; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgy = &gy[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                RmgType *bgz = &gz[(iy - order/2)*dimz + (ix - order/2)*dimy*dimz - order/2];
                for (int iz = order/2; iz < dimz + order/2; iz++)
                {
                    bgz[iz] += -cz[0] * A[iz + 1*ixs - 1*iys - 1] + cz[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        bgz[iz] += -cz[1] * A[iz + 2*ixs - 2*iys - 2] + cz[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        bgz[iz] += -cz[2] * A[iz + 3*ixs - 3*iys - 3] + cz[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        bgz[iz] += -cz[3] * A[iz + 4*ixs - 4*iys - 4] + cz[3] * A[iz - 4*ixs + 4*iys + 4];
                    if constexpr(order >= 10)
                        bgz[iz] += -cz[4] * A[iz + 5*ixs - 5*iys - 5] + cz[4] * A[iz - 5*ixs + 5*iys + 5];
                    if constexpr(order >= 12)
                        bgz[iz] += -cz[5] * A[iz + 6*ixs - 6*iys - 6] + cz[5] * A[iz - 6*ixs + 6*iys + 6];

                    bgy[iz] += -cy[0] * A[iz + 1*ixs - 1*iys - 1] + cy[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs - 2*iys - 2] + cy[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs - 3*iys - 3] + cy[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs - 4*iys - 4] + cy[3] * A[iz - 4*ixs + 4*iys + 4];
                    if constexpr(order >= 10)
                        bgy[iz] += -cy[4] * A[iz + 5*ixs - 5*iys - 5] + cy[4] * A[iz - 5*ixs + 5*iys + 5];
                    if constexpr(order >= 12)
                        bgy[iz] += -cy[5] * A[iz + 6*ixs - 6*iys - 6] + cy[5] * A[iz - 6*ixs + 6*iys + 6];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs - 1*iys - 1] + cx[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs - 2*iys - 2] + cx[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs - 3*iys - 3] + cx[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs - 4*iys - 4] + cx[3] * A[iz - 4*ixs + 4*iys + 4];
                    if constexpr(order >= 10)
                        bgx[iz] += -cx[4] * A[iz + 5*ixs - 5*iys - 5] + cx[4] * A[iz - 5*ixs + 5*iys + 5];
                    if constexpr(order >= 12)
                        bgx[iz] += -cx[5] * A[iz + 6*ixs - 6*iys - 6] + cx[5] * A[iz - 6*ixs + 6*iys + 6];
                }                   /* end for */
            }
        }
    }
} /* end app8_gradient_general */



// Gets the central coefficient
double FiniteDiff::fd_coeff0(int order, double hxgrid)
{
    LaplacianCoeff *LC1, *LC2;
    LC2 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order];
    if(order == 2)
        LC1 = LC2;
    else
        LC1 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order-2];

    double scale = LC2->gen_hxgrid / hxgrid;
    scale = scale*scale;

    double coeff0 = 0.0;
    for(int ax=0;ax < 13;ax++)
    {
        double c1, c2 = 0.0;
        if(this->alt_laplacian) c2 = cfac[0];
        if(order == 2) {c1=1.0;c2 = 0.0;}    // no optimzation for 2nd order
        c1 = 1.0 + c2;
        coeff0 += c1*LC2->plane_centers[ax] - c2*LC1->plane_centers[ax];
    }
    return scale*coeff0;
}


// Computes combined coefficients
template <typename RmgType>
void FiniteDiff::fd_combined_coeffs(int order, double hxgrid, int ax, RmgType * cm, RmgType *cp, double *kvec)
{
    double s1 = 2.0;
    RmgType t1;
    RmgType x1, y1, z1;
    LaplacianCoeff *LC1, *LC2;
    LC2 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order];
    if(order == 2)
        LC1 = LC2;
    else
        LC1 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order-2];

    double scale = LC2->gen_hxgrid / hxgrid;
    scale = scale*scale;

    RmgType I_t;
    if(typeid(RmgType) == typeid(double))
    {
        I_t = 0.0;
    }
    else if(typeid(RmgType) == typeid(float))
    {
        I_t = 0.0;
    }
    else if(typeid(RmgType) == typeid(std::complex<double>))
    {
        double *iptr = (double *)&I_t;
        iptr[0] = 0.0;iptr[1] = 1.0;
    }
    else if(typeid(RmgType) == typeid(std::complex<float>))
    {
        float *iptr = (float *)&I_t;
        iptr[0] = 0.0;iptr[1] = 1.0;
    }

    // 2nd order if for multigrid and has no optimizations
    if(order == 2)
    {
        t1 = scale*LC2->axis_lc[ax][0];
        x1 = scale*LC2->axis_gc_x[ax][0];
        y1 = scale*LC2->axis_gc_y[ax][0];
        z1 = scale*LC2->axis_gc_z[ax][0];
        cm[0] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
        cp[0] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
        return;
    }

    double c1, c2=0.0, c1g, c2g=0.0;
    if(this->alt_laplacian && order > 2) c2 = cfac[0];
    if(this->alt_laplacian && order > 2) c2g = cfac[1];
    c1 = scale*(1.0 + c2);
    c1g = scale*(1.0 + c2g);
    //t1 = c1*LC2->axis_lc[ax][3] - c2*LC1->axis_lc[ax][2];

    for(int i=0;i < order/2-1;i++)
    {
        t1 = c1*LC2->axis_lc[ax][order/2-i-1] - c2*LC1->axis_lc[ax][order/2-i-2];
        x1 = c1g*LC2->axis_gc_x[ax][order/2-i-1] - c2g*LC1->axis_gc_x[ax][order/2-i-2];
        y1 = c1g*LC2->axis_gc_y[ax][order/2-i-1] - c2g*LC1->axis_gc_y[ax][order/2-i-2];
        z1 = c1g*LC2->axis_gc_z[ax][order/2-i-1] - c2g*LC1->axis_gc_z[ax][order/2-i-2];
        cm[i] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
        cp[i] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    }
    t1 = c1*LC2->axis_lc[ax][0];
    x1 = c1g*LC2->axis_gc_x[ax][0];
    y1 = c1g*LC2->axis_gc_y[ax][0];
    z1 = c1g*LC2->axis_gc_z[ax][0];
    cm[order/2-1] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[order/2-1] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

}

template <typename RmgType>
void FiniteDiff::fd_gradient_coeffs(int order, double hxgrid, int axis , RmgType *cx, RmgType *cy, RmgType *cz)
{
    LaplacianCoeff *LC1, *LC2;
    LC2 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order];
    if(order == 2)
        LC1 = LC2;
    else
        LC1 = FiniteDiff::FdCoeffs[LCkey(hxgrid)+order-2];

    double scale = LC2->gen_hxgrid / hxgrid;
    scale = scale*scale;

    double c1, c2=0.0;
    if(this->alt_laplacian && order > 2) c2 = cfac[1];
    c1 = scale*(1.0 + c2);

    for(int i=0;i < order/2-1;i++)
    {
        cx[i] = c1*LC2->axis_gc_x[axis][order/2-i-1] - c2*LC1->axis_gc_x[axis][order/2-i-2];
        cy[i] = c1*LC2->axis_gc_y[axis][order/2-i-1] - c2*LC1->axis_gc_y[axis][order/2-i-2];
        cz[i] = c1*LC2->axis_gc_z[axis][order/2-i-1] - c2*LC1->axis_gc_z[axis][order/2-i-2];
    }
    cx[order/2-1] = c1*LC2->axis_gc_x[axis][0];
    cy[order/2-1] = c1*LC2->axis_gc_y[axis][0];
    cz[order/2-1] = c1*LC2->axis_gc_z[axis][0];
}


