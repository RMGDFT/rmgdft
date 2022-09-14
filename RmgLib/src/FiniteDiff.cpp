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
#include "LaplacianCoeff.h"
#include "RmgTimer.h"
#include "rmg_error.h"

#define         PI          3.14159265358979323

int FiniteDiff::allocation_limit = 65536;
double FiniteDiff::cfac[12];

// Force instantiation of float, double and complex versions.
template double FiniteDiff::app_cil_sixth<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_sixth<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_sixth<std::complex<float> >(std::complex <float> *, std::complex <float> *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_sixth<std::complex<double> >(std::complex <double> *, std::complex <double> *, int, int, int, double, double, double);

template double FiniteDiff::app_del2_np<double>(double *, double *, double, double, double);
template void FiniteDiff::app_cir_sixth<float>(float *, float *, int, int, int);
template void FiniteDiff::app_cir_sixth<double>(double *, double *, int, int, int);
template void FiniteDiff::app_cir_sixth<std::complex<float> >(std::complex <float> *, std::complex <float> *, int, int, int);
template void FiniteDiff::app_cir_sixth<std::complex<double> >(std::complex <double> *, std::complex <double> *, int, int, int);


template double FiniteDiff::app_cil_fourth<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth<std::complex<float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth<std::complex<double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template double FiniteDiff::app_cil_fourth_threaded<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth_threaded<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth_threaded<std::complex<float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app_cil_fourth_threaded<std::complex<double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template void FiniteDiff::app_cir_fourth<float>(float *, float *, int, int, int);
template void FiniteDiff::app_cir_fourth<double>(double *, double *, int, int, int);
template void FiniteDiff::app_cir_fourth<std::complex<float> >(std::complex<float> *, std::complex<float> *, int, int, int);
template void FiniteDiff::app_cir_fourth<std::complex<double> >(std::complex<double> *, std::complex<double> *, int, int, int);

template double FiniteDiff::app2_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app2_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template double FiniteDiff::app2_del2_offset<float>(float *, float *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<double>(double *, double *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double FiniteDiff::app2_del2_offset<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template double FiniteDiff::app8_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app8_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template double FiniteDiff::app10_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template void FiniteDiff::app_gradient_sixth<float> (float *, float *, float *, float *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_sixth<double> (double *, double *, double *, double *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_sixth<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int, double , double , double );
template void FiniteDiff::app_gradient_sixth<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int, double , double , double );

template void FiniteDiff::app_gradient_eighth<float> (float *, float *, float *, float *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_eighth<double> (double *, double *, double *, double *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_eighth<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int, double , double , double );
template void FiniteDiff::app_gradient_eighth<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int, double , double , double );

template void FiniteDiff::app_gradient_tenth<float> (float *, float *, float *, float *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_tenth<double> (double *, double *, double *, double *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_tenth<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int, double , double , double );
template void FiniteDiff::app_gradient_tenth<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int, double , double , double );

template double FiniteDiff::app8_combined<float>(float *, float *, int, int, int, double, double, double, double *kvec);
template double FiniteDiff::app8_combined<double>(double *, double *, int, int, int, double, double, double, double *kvec);
template double FiniteDiff::app8_combined<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec);
template double FiniteDiff::app8_combined<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec);
template double FiniteDiff::app8_combined<float>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<double>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

template double FiniteDiff::app8_combined_orthorhombic<float>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_orthorhombic<double>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_orthorhombic<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_orthorhombic<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

template double FiniteDiff::app8_combined_general<float>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_general<double>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_general<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined_general<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

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
double FiniteDiff::app_cil_sixth (RmgType *rptr, RmgType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    RmgTimer RT("App_cil: computation");
    int iz, ix, iy, incx, incy, incxr, incyr, ibrav;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    double ihx, ihy, ihz;
    RmgType rz, rzms, rzps, rzpps;
    RmgType rfc1, rbc1, rbc2, rd1, rd2, rd3, rd4;
    RmgType td1, td2, td3, td4, td5, td6, td7, td8, tdx;
    double c0 = -116.0 / 90.0;
    double c1 = 31.0 / 232.0;
    double c2 = 49.0 / 60.0;
    double c3 = -31.0 / 464.0;
    double c4 = 1.0 / 10.0;
    double c5 = 1.0 / 120.0;
    double c6 = -1.0 / 240.0;
    double c7 = 1.0 / 144.0;

    ibrav = L->get_ibrav_type();

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    ihx = 1.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    ihy = 1.0 / (gridhy * gridhy * L->get_yside() * L->get_yside());
    ihz = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());

    double ccd (c0 * (ihx + ihy + ihz));
    RmgType cc (std::real(c0 * (ihx + ihy + ihz)));
    RmgType fcx ( std::real(c1 * ccd + c2 * ihx));
    RmgType fcy ( std::real(c1 * ccd + c2 * ihy));
    RmgType fcz ( std::real(c1 * ccd + c2 * ihz));

    RmgType ecxy ( std::real(c3 * ccd - c4 * ihz));
    RmgType ecxz ( std::real(c3 * ccd - c4 * ihy));
    RmgType ecyz ( std::real(c3 * ccd - c4 * ihx));

    RmgType cor ( std::real(c7 * (ihx + ihy + ihz)));

    RmgType fc2x ( std::real(c5 * (ihy + ihz)));
    RmgType fc2y ( std::real(c5 * (ihx + ihz)));
    RmgType fc2z ( std::real(c5 * (ihx + ihy)));

    RmgType tcx ( std::real(c6 * ihx));
    RmgType tcy ( std::real(c6 * ihy));
    RmgType tcz ( std::real(c6 * ihz));

    // Handle the general case first
    if((dimz % 4) || (ibrav != CUBIC_PRIMITIVE)) {

        for (ix = 2; ix < dimx + 2; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            ixmms = (ix - 2) * incx;
            ixpps = (ix + 2) * incx;

            for (iy = 2; iy < dimy + 2; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                iymms = (iy - 2) * incy;
                iypps = (iy + 2) * incy;

                for (iz = 2; iz < dimz + 2; iz++)
                {

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                        fcx * (rptr[ixms + iys + iz] + rptr[ixps + iys + iz]) +
                        fcy * (rptr[ixs + iyms + iz] + rptr[ixs + iyps + iz]) +
                        fcz * (rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        ecxz * (rptr[ixms + iys + (iz - 1)] + rptr[ixps + iys + (iz - 1)] +
                                rptr[ixms + iys + (iz + 1)] + rptr[ixps + iys + (iz + 1)]) +
                        ecyz * (rptr[ixs + iyms + (iz - 1)] + rptr[ixs + iyps + (iz - 1)] +
                                rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyps + (iz + 1)]) +
                        ecxy * (rptr[ixms + iyms + iz] + rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] + rptr[ixps + iyps + iz]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        cor * (rptr[ixms + iyms + (iz - 1)] + rptr[ixps + iyms + (iz - 1)] +
                               rptr[ixms + iyps + (iz - 1)] + rptr[ixps + iyps + (iz - 1)] +
                               rptr[ixms + iyms + (iz + 1)] + rptr[ixps + iyms + (iz + 1)] +
                               rptr[ixms + iyps + (iz + 1)] + rptr[ixps + iyps + (iz + 1)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz]) +
                        fc2y * (rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                        fc2z * (rptr[ixs + iys + (iz - 2)] + rptr[ixs + iys + (iz + 2)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        tcx * (rptr[ixps + iypps + iz] + rptr[ixps + iymms + iz] +
                               rptr[ixms + iypps + iz] + rptr[ixms + iymms + iz] +
                               rptr[ixps + iys + (iz + 2)] + rptr[ixps + iys + (iz - 2)] +
                               rptr[ixms + iys + (iz + 2)] + rptr[ixms + iys + (iz - 2)]) +
                        tcy * (rptr[ixpps + iyps + iz] + rptr[ixmms + iyps + iz] +
                               rptr[ixpps + iyms + iz] + rptr[ixmms + iyms + iz] +
                               rptr[ixs + iyps + (iz + 2)] + rptr[ixs + iyps + (iz - 2)] +
                               rptr[ixs + iyms + (iz + 2)] + rptr[ixs + iyms + (iz - 2)]) +
                        tcz * (rptr[ixpps + iys + (iz + 1)] + rptr[ixmms + iys + (iz + 1)] +
                               rptr[ixpps + iys + (iz - 1)] + rptr[ixmms + iys + (iz - 1)] +
                               rptr[ixs + iypps + (iz + 1)] + rptr[ixs + iymms + (iz + 1)] +
                               rptr[ixs + iypps + (iz - 1)] + rptr[ixs + iymms + (iz - 1)]);

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

        return (double)std::real(cc);

    }

    // Optimized case for dimz divisible by 4 and cubic primitive grid

    for (ix = 2; ix < dimx + 2; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < dimy + 2; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            // Compute the middle set of edges (2nd nn) before entry to the loop
            rz = rptr[ixms + iys + 2] +
                 rptr[ixps + iys + 2] +
                 rptr[ixs + iyms + 2] +
                 rptr[ixs + iyps + 2];

            // Compute the trailing set of edges (2nd nn) before entry to the  loop
            rzms = rptr[ixps + iys + 1] +
                   rptr[ixms + iys + 1] +
                   rptr[ixs + iyps + 1] +
                   rptr[ixs + iyms + 1];

            // Compute the 4 nn with same z value as loop index on entry
            rfc1 = rptr[ixms + iyms + 2] + rptr[ixms + iyps + 2] +
                   rptr[ixps + iyms + 2] + rptr[ixps + iyps + 2];

            // Compute a pair trailing sets of corners on entry
            rbc1 = rptr[ixms + iyms + 1] + rptr[ixps + iyms + 1] +
                   rptr[ixms + iyps + 1] + rptr[ixps + iyps + 1];
            rbc2 = rptr[ixms + iyms + 2] + rptr[ixps + iyms + 2] +
                   rptr[ixms + iyps + 2] + rptr[ixps + iyps + 2];

            rd1 = rptr[ixpps + iys + 1] + rptr[ixmms + iys + 1] +
                  rptr[ixs + iypps + 1] + rptr[ixs + iymms + 1];

            rd4 = rptr[ixpps + iys + 2] + rptr[ixmms + iys + 2] +
                  rptr[ixs + iypps + 2] + rptr[ixs + iymms + 2];


            td2 =   rptr[ixps + iys] +
                           rptr[ixms + iys] +
                           rptr[ixs + iyps] +
                           rptr[ixs + iyms];

            td4 =   rptr[ixps + iys + 1] +
                           rptr[ixms + iys + 1] +
                           rptr[ixs + iyps + 1] +
                           rptr[ixs + iyms + 1];

            td6 =   rptr[ixps + iys + 2] +
                           rptr[ixms + iys + 2] +
                           rptr[ixs + iyps + 2] +
                           rptr[ixs + iyms + 2];

            td8 =   rptr[ixps + iys + 3] +
                           rptr[ixms + iys + 3] +
                           rptr[ixs + iyps + 3] +
                           rptr[ixs + iyms + 3];

            td1 =   rptr[ixps + iys + 4] +
                           rptr[ixms + iys + 4] +
                           rptr[ixs + iyps + 4] +
                           rptr[ixs + iyms + 4];

            td3 =   rptr[ixps + iys + 5] +
                           rptr[ixms + iys + 5] +
                           rptr[ixs + iyps + 5] +
                           rptr[ixs + iyms + 5];


            td5 =   rptr[ixps + iys + 6] +
                           rptr[ixms + iys + 6] +
                           rptr[ixs + iyps + 6] +
                           rptr[ixs + iyms + 6];


            td7 =   rptr[ixps + iys + 7] +
                           rptr[ixms + iys + 7] +
                           rptr[ixs + iyps + 7] +
                           rptr[ixs + iyms + 7];



            for (iz = 2; iz < dimz + 2; iz+=4)
            {

                tdx =   rptr[ixps + iys + iz + 6] +
                        rptr[ixms + iys + iz + 6] +
                        rptr[ixs + iyps + iz + 6] +
                        rptr[ixs + iyms + iz + 6];

                // Forward set of edges
                rzps = rptr[ixps + iys + iz + 1] +
                       rptr[ixms + iys + iz + 1] +
                       rptr[ixs + iyps + iz + 1] +
                       rptr[ixs + iyms + iz + 1];

                // First
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                    fcx * rz +
                    fcx * (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 1]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    ecxz * rzms + ecxz * rzps + ecxz * rfc1;

                // Compute the forward set of corners
                rfc1 = rptr[ixms + iyms + iz + 1] + rptr[ixps + iyms + iz + 1] +
                       rptr[ixms + iyps + iz + 1] + rptr[ixps + iyps + iz + 1];

                rd3 = rptr[ixpps + iys + iz + 1] + rptr[ixmms + iys + iz + 1] +
                          rptr[ixs + iypps + iz + 1] + rptr[ixs + iymms + iz + 1];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    cor * rbc1 + cor * rfc1;
                rbc1 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                    fc2x * (rptr[ixs + iys + iz - 2] + rptr[ixs + iys + iz + 2]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    tcx * rptr[ixps + iypps + iz] + tcx * rptr[ixps + iymms + iz] +
                    tcx * rptr[ixms + iypps + iz] + tcx * rptr[ixms + iymms + iz] +
                    tcx * rptr[ixpps + iyps + iz] + tcx * rptr[ixmms + iyps + iz] +
                    tcx * rptr[ixpps + iyms + iz] + tcx * rptr[ixmms + iyms + iz] +
                    tcx * rd1 +
                    tcx * rd3 +
                    tcx * (td1 + td2);

                td2 = td1;
                td1 = tdx;
                tdx =   rptr[ixps + iys + iz + 7] +
                        rptr[ixms + iys + iz + 7] +
                        rptr[ixs + iyps + iz + 7] +
                        rptr[ixs + iyms + iz + 7];

                // Compute another set of forward edges here
                rzpps = rptr[ixps + iys + iz + 2] +
                        rptr[ixms + iys + iz + 2] +
                        rptr[ixs + iyps + iz + 2] +
                        rptr[ixs + iyms + iz + 2];

                // Second
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] = cc * rptr[ixs + iys + iz + 1] +
                    fcx * rzps +
                    fcx * (rptr[ixs + iys + iz] + rptr[ixs + iys + iz + 2]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    ecxz * rz +  ecxz * rzpps +  ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 2] + rptr[ixps + iyms + iz + 2] +
                       rptr[ixms + iyps + iz + 2] + rptr[ixps + iyps + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    cor * rbc2 + cor * rfc1;
                rbc2 = rfc1;

                rd2 = rptr[ixpps + iys + iz + 2] + rptr[ixmms + iys + iz + 2] +
                      rptr[ixs + iypps + iz + 2] + rptr[ixs + iymms + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    fc2x * (rptr[ixmms + iys + iz + 1] + rptr[ixpps + iys + iz + 1] +
                            rptr[ixs + iymms + iz + 1] + rptr[ixs + iypps + iz + 1]) +
                    fc2x *  (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 3]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    tcx * rptr[ixps + iypps + iz + 1] + tcx * rptr[ixps + iymms + iz + 1] +
                    tcx * rptr[ixms + iypps + iz + 1] + tcx * rptr[ixms + iymms + iz + 1] +
                    tcx * rptr[ixpps + iyps + iz + 1] + tcx * rptr[ixmms + iyps + iz + 1] +
                    tcx * rptr[ixpps + iyms + iz + 1] + tcx * rptr[ixmms + iyms + iz + 1] +
                    tcx * rd4 +
                    tcx * rd2 +
                    tcx * (td3 + td4);

                td4 = td3;
                td3 = tdx;
                tdx =   rptr[ixps + iys + iz + 8] +
                           rptr[ixms + iys + iz + 8] +
                           rptr[ixs + iyps + iz + 8] +
                           rptr[ixs + iyms + iz + 8];

                // Compute the forward set of edges here which becomes the trailing set at the top
                rzms = rptr[ixms + iys + iz+3] +
                       rptr[ixps + iys + iz+3] +
                       rptr[ixs + iyms + iz+3] +
                       rptr[ixs + iyps + iz+3];

                // Third
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] = cc * rptr[ixs + iys + iz + 2] +
                    fcx * rzpps +
                    fcx * (rptr[ixs + iys + iz + 1] + rptr[ixs + iys + iz + 3]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    ecxz * rzps + ecxz * rzms + ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 3] + rptr[ixps + iyms + iz + 3] +
                       rptr[ixms + iyps + iz + 3] + rptr[ixps + iyps + iz + 3];

                rd1 = rptr[ixpps + iys + iz + 3] + rptr[ixmms + iys + iz + 3] +
                      rptr[ixs + iypps + iz + 3] + rptr[ixs + iymms + iz + 3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    cor * rbc1 + cor * rfc1;
                rbc1 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    fc2x * (rptr[ixmms + iys + iz + 2] + rptr[ixpps + iys + iz + 2] +
                            rptr[ixs + iymms + iz + 2] + rptr[ixs + iypps + iz + 2]) +
                    fc2x *  (rptr[ixs + iys + iz] + rptr[ixs + iys + iz + 4]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz)] +=
                    tcx * rptr[ixps + iypps + iz + 2] + tcx * rptr[ixps + iymms + iz + 2] +
                    tcx * rptr[ixms + iypps + iz + 2] + tcx * rptr[ixms + iymms + iz + 2] +
                    tcx * rptr[ixpps + iyps + iz + 2] + tcx * rptr[ixmms + iyps + iz + 2] +
                    tcx * rptr[ixpps + iyms + iz + 2] + tcx * rptr[ixmms + iyms + iz + 2] +
                    tcx * rd3 +
                    tcx * rd1 +
                    tcx * (td5 + td6);

                td6 = td5;
                td5 = tdx;
                tdx =   rptr[ixps + iys + iz + 9] +
                        rptr[ixms + iys + iz + 9] +
                        rptr[ixs + iyps + iz + 9] +
                        rptr[ixs + iyms + iz + 9];

                // Compute the forward set of edges here which becomes the middle set at the top
                rz = rptr[ixps + iys + iz + 4] +
                     rptr[ixms + iys + iz + 4] +
                     rptr[ixs + iyps + iz + 4] +
                     rptr[ixs + iyms + iz + 4];

                // Fourth
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] = cc * rptr[ixs + iys + iz + 3] +
                    fcx * rzms +
                    fcx * (rptr[ixs + iys + iz + 2] + rptr[ixs + iys + iz + 4]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    ecxz * rzpps + ecxz * rz + ecxz * rfc1;

                // Another set of forward corners here
                rfc1 = rptr[ixms + iyms + iz + 4] + rptr[ixps + iyms + iz + 4] +
                       rptr[ixms + iyps + iz + 4] + rptr[ixps + iyps + iz + 4];

                rd4 = rptr[ixpps + iys + iz + 4] + rptr[ixmms + iys + iz + 4] +
                      rptr[ixs + iypps + iz + 4] + rptr[ixs + iymms + iz + 4];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    cor * rbc2 + cor * rfc1;
                rbc2 = rfc1;

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    fc2x * (rptr[ixmms + iys + iz + 3] + rptr[ixpps + iys + iz + 3] +
                            rptr[ixs + iymms + iz + 3] + rptr[ixs + iypps + iz + 3]) +
                    fc2x * (rptr[ixs + iys + iz + 1] + rptr[ixs + iys + iz + 5]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz + 1)] +=
                    tcx * rptr[ixps + iypps + iz + 3] + tcx * rptr[ixps + iymms + iz + 3] +
                    tcx * rptr[ixms + iypps + iz + 3] + tcx * rptr[ixms + iymms + iz + 3] +
                    tcx * rptr[ixpps + iyps + iz + 3] + tcx * rptr[ixmms + iyps + iz + 3] +
                    tcx * rptr[ixpps + iyms + iz + 3] + tcx * rptr[ixmms + iyms + iz + 3] +
                    tcx * rd2 +
                    tcx * rd4 +
                    tcx * (td7 + td8);

                td8 = td7;
                td7 = tdx;



            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    return (double)std::real(cc);

}



template <typename RmgType>
void FiniteDiff::app_cir_sixth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz)
{

    RmgTimer RT("App_cir: computation");
    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx, ixmms, ixpps, iymms, iypps;
    int incyr, incxr;
    RmgType rz, rzps, rzms, rzpps;
    RmgType c000, c100, c110, c200;

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;

    // Handle the general case first
    if(dimz % 4) {

        for (ix = 2; ix < dimx + 2; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            ixmms = (ix - 2) * incx;
            ixpps = (ix + 2) * incx;

            for (iy = 2; iy < dimy + 2; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                iymms = (iy - 2) * incy;
                iypps = (iy + 2) * incy;

                for (iz = 2; iz < dimz + 2; iz++)
                {

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                        c100 * (rptr[ixs + iys + (iz - 1)] +
                                rptr[ixs + iys + (iz + 1)] +
                                rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        c110 * (rptr[ixps + iyps + iz] +
                                rptr[ixps + iyms + iz] +
                                rptr[ixms + iyps + iz] +
                                rptr[ixms + iyms + iz] +
                                rptr[ixps + iys + (iz + 1)] +
                                rptr[ixps + iys + (iz - 1)] +
                                rptr[ixms + iys + (iz + 1)] +
                                rptr[ixms + iys + (iz - 1)] +
                                rptr[ixs + iyps + (iz + 1)] +
                                rptr[ixs + iyps + (iz - 1)] +
                                rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyms + (iz - 1)]);

                    b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                        c200 * (rptr[ixs + iys + (iz - 2)] +
                                rptr[ixs + iys + (iz + 2)] +
                                rptr[ixmms + iys + iz] +
                                rptr[ixpps + iys + iz] +
                                rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]);

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

        return;
    }


    // Optimized case for dimz divisible by 4
    for (ix = 2; ix < dimx + 2; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < dimy + 2; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            // Compute the middle set of edges before entry to the loop
            rz = rptr[ixms + iys + 2] +
                 rptr[ixps + iys + 2] +
                 rptr[ixs + iyms + 2] +
                 rptr[ixs + iyps + 2];

            // Compute the trailing set of edges before entry to the  loop
            rzms = rptr[ixps + iys + 1] +
                   rptr[ixms + iys + 1] +
                   rptr[ixs + iyps + 1] +
                   rptr[ixs + iyms + 1];


            for (iz = 2; iz < dimz + 2; iz+=4)
            {

                // Forward set of edges
                rzps = rptr[ixps + iys + iz + 1] +
                       rptr[ixms + iys + iz + 1] +
                       rptr[ixs + iyps + iz + 1] +
                       rptr[ixs + iyms + iz + 1];

                // First
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (rptr[ixs + iys + iz - 1] +
                            rptr[ixs + iys + iz + 1] +
                            rz) + c000 * rptr[ixs + iys + iz];


                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c110 * (rptr[ixps + iyps + iz] +
                            rptr[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz] +
                            rzms +
                            rzps);

                // Compute another set of forward edges here
                rzpps = rptr[ixps + iys + iz + 2] +
                        rptr[ixms + iys + iz + 2] +
                        rptr[ixs + iyps + iz + 2] +
                        rptr[ixs + iyms + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (rptr[ixs + iys + iz - 2] +
                            rptr[ixs + iys + iz + 2]) +
                    c200 * (rptr[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + 
                            rptr[ixs + iypps + iz]);


                // Second
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz)] +
                            rptr[ixs + iys + iz + 2] +
                            rzps) + 
                            c000 * rptr[ixs + iys + iz + 1];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    c110 * (rptr[ixps + iyps + iz+1] +
                            rptr[ixps + iyms + iz+1] +
                            rptr[ixms + iyps + iz+1] +
                            rptr[ixms + iyms + iz+1] +
                            rzpps +
                            rz);

                // Compute the forward set of edges here which becomes the trailing set at the top
                rzms = rptr[ixms + iys + iz+3] +
                       rptr[ixps + iys + iz+3] +
                       rptr[ixs + iyms + iz+3] +
                       rptr[ixs + iyps + iz+3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 1)] +=
                    c200 * (rptr[ixs + iys + iz - 1] +
                            rptr[ixs + iys + iz + 3]) +
                    c200 * (rptr[ixmms + iys + iz+1] +
                            rptr[ixpps + iys + iz+1] +
                            rptr[ixs + iymms + iz+1] + 
                            rptr[ixs + iypps + iz+1]);

                // Third
                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] =
                    c100 * (rptr[ixs + iys + iz+1] +
                            rptr[ixs + iys + iz + 3] +
                            rzpps) +
                            c000 * rptr[ixs + iys + iz + 2];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    c110 * (rptr[ixps + iyps + iz+2] +
                            rptr[ixps + iyms + iz+2] +
                            rptr[ixms + iyps + iz+2] +
                            rptr[ixms + iyms + iz+2] +
                            rzms +
                            rzps);

                // Compute the forward set of edges here which becomes the middle set at the top
                rz = rptr[ixps + iys + iz + 4] +
                     rptr[ixms + iys + iz + 4] +
                     rptr[ixs + iyps + iz + 4] +
                     rptr[ixs + iyms + iz + 4];

                b[(ix - 2) * incxr + (iy - 2) * incyr + iz] +=
                    c200 * (rptr[ixs + iys + (iz)] +
                            rptr[ixs + iys + (iz + 4)]) +
                    c200 * (rptr[ixmms + iys + iz+2] +
                            rptr[ixpps + iys + iz+2] +
                            rptr[ixs + iymms + iz+2] + 
                            rptr[ixs + iypps + iz+2]);

                // Fourth
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] =
                    c100 * (rptr[ixs + iys + iz+2] +
                            rptr[ixs + iys + iz + 4] +
                            rzms) + 
                            c000 * rptr[ixs + iys + iz + 3];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] +=
                    c110 * (rptr[ixps + iyps + iz+3] +
                            rptr[ixps + iyms + iz+3] +
                            rptr[ixms + iyps + iz+3] +
                            rptr[ixms + iyms + iz+3] +
                            rzpps +
                            rz);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz+1)] +=
                    c200 * (rptr[ixs + iys + iz+1] +
                            rptr[ixs + iys + iz + 5]) +
                    c200 * (rptr[ixmms + iys + iz+3] +
                            rptr[ixpps + iys + iz+3] +
                            rptr[ixs + iymms + iz+3] + 
                            rptr[ixs + iypps + iz+3]);


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

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



template <typename RmgType>
double FiniteDiff::app8_del2(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
               double gridhx, double gridhy, double gridhz)
{
    double kvec[3] = {0.0,0.0,0.0};
    return FiniteDiff::app8_combined_general(a, b, dimx, dimy, dimz,
        gridhx, gridhy, gridhz, kvec, false);

}  /* end app8_del2 */


template <typename RmgType>
double FiniteDiff::app10_del2(RmgType * a, RmgType * b, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{

    int ibrav = L->get_ibrav_type();
    int iz, ix, iy;
    int ixs, iys, ix1, iy1;

    ixs = (dimy + 10) * (dimz + 10);
    iys = (dimz + 10);
    ix1 = dimy * dimz;
    iy1 = dimz;

    // nine and seven point stencils, 2nd derivative, extrapolated
    int ic = 5;
    double x[12], w1[12], w2[12];
    for(int i=0;i<12;i++) x[i] = (double)i;
    gen_weights(11, 2, (double)ic, x, w1);
    gen_weights(9, 2, (double)(ic-1), x, w2);
    double hf = 1.0, c1, c2=0.0;
    if(ibrav == HEXAGONAL || ibrav == HEXAGONAL2) hf = 2.0/3.0;

    double d1 = -8.0/3150.0;
    double d2 = 6.01250601250605e-5;
    double dr = d1 / d2;
    double k2 = PI*PI/8.0;
    double h2x = gridhx * gridhx * L->get_xside() * L->get_xside();
    double h2y = gridhy * gridhy * L->get_yside() * L->get_yside();
    double h2z = gridhz * gridhz * L->get_zside() * L->get_zside();
    double maxh = std::max(h2x, h2y);
    maxh = std::max(maxh, h2z);

    double hadj = sqrt(h2x / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    double th2 = (c1*w1[ic] - c2*w2[ic-1]) / h2x;
    RmgType t1x ((c1*w1[ic+1] - c2*w2[ic]) * hf / h2x);
    RmgType t2x ((c1*w1[ic+2] - c2*w2[ic+1]) * hf / h2x);
    RmgType t3x ((c1*w1[ic+3] - c2*w2[ic+2]) * hf / h2x);
    RmgType t4x ((c1*w1[ic+4] - c2*w2[ic+3]) * hf / h2x);
    RmgType t5x (c1*w1[ic+5] * hf / h2x);

    hadj = sqrt(h2y / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) /  h2y;
    RmgType t1y ((c1*w1[ic+1] - c2*w2[ic]) / h2y);
    RmgType t2y ((c1*w1[ic+2] - c2*w2[ic+1]) / h2y);
    RmgType t3y ((c1*w1[ic+3] - c2*w2[ic+2]) / h2y);
    RmgType t4y ((c1*w1[ic+4] - c2*w2[ic+3]) / h2y);
    RmgType t5y (c1*w1[ic+5] / h2y);

    hadj = sqrt(h2z / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hf*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) /  h2z;
    RmgType t1z ((c1*w1[ic+1] - c2*w2[ic]) / h2z);
    RmgType t2z ((c1*w1[ic+2] - c2*w2[ic+1]) / h2z);
    RmgType t3z ((c1*w1[ic+3] - c2*w2[ic+2]) / h2z);
    RmgType t4z ((c1*w1[ic+4] - c2*w2[ic+3]) / h2z);
    RmgType t5z (c1*w1[ic+5] / h2z);
    RmgType t0 (th2);

    // NULL b means we just want the diagonal component.
    if(b == NULL) return (double)std::real(t0);

    switch(ibrav)
    {
        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (ix = 5; ix < dimx + 5; ix++)
            {

                for (iy = 5; iy < dimy + 5; iy++)
                {

                    for (iz = 5; iz < dimz + 5; iz++)
                    {

                        b[(ix - 5) * ix1 + (iy - 5) * iy1 + iz - 5] =
                            t0 * a[ix * ixs + iy * iys + iz] +
                            t1x * a[(ix - 1) * ixs + iy * iys + iz] +
                            t1x * a[(ix + 1) * ixs + iy * iys + iz] +
                            t2x * a[(ix - 2) * ixs + iy * iys + iz] +
                            t2x * a[(ix + 2) * ixs + iy * iys + iz] +
                            t3x * a[(ix - 3) * ixs + iy * iys + iz] +
                            t3x * a[(ix + 3) * ixs + iy * iys + iz] +
                            t4x * a[(ix - 4) * ixs + iy * iys + iz] +
                            t4x * a[(ix + 4) * ixs + iy * iys + iz] +
                            t5x * a[(ix - 5) * ixs + iy * iys + iz] +
                            t5x * a[(ix + 5) * ixs + iy * iys + iz] +

                            t1y * a[ix * ixs + (iy - 1) * iys + iz] +
                            t1y * a[ix * ixs + (iy + 1) * iys + iz] +
                            t2y * a[ix * ixs + (iy - 2) * iys + iz] +
                            t2y * a[ix * ixs + (iy + 2) * iys + iz] +
                            t3y * a[ix * ixs + (iy - 3) * iys + iz] +
                            t3y * a[ix * ixs + (iy + 3) * iys + iz] +
                            t4y * a[ix * ixs + (iy - 4) * iys + iz] +
                            t4y * a[ix * ixs + (iy + 4) * iys + iz] +
                            t5y * a[ix * ixs + (iy - 5) * iys + iz] +
                            t5y * a[ix * ixs + (iy + 5) * iys + iz] +

                            t1z * a[ix * ixs + iy * iys + iz - 1] +
                            t1z * a[ix * ixs + iy * iys + iz + 1] +
                            t2z * a[ix * ixs + iy * iys + iz - 2] +
                            t2z * a[ix * ixs + iy * iys + iz + 2] +
                            t3z * a[ix * ixs + iy * iys + iz - 3] +
                            t3z * a[ix * ixs + iy * iys + iz + 3] +
                            t4z * a[ix * ixs + iy * iys + iz - 4] +
                            t4z * a[ix * ixs + iy * iys + iz + 4] +
                            t5z * a[ix * ixs + iy * iys + iz - 5] +
                            t5z * a[ix * ixs + iy * iys + iz + 5];

                    }                   /* end for */
                }                       /* end for */
            }                           /* end for */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            for (int ix = 5; ix < dimx + 5; ix++)
            {

                for (int iy = 5; iy < dimy + 5; iy++)
                {

                    RmgType *A = &a[iy*iys + ix*ixs];
                    RmgType *B = &b[(iy - 5)*dimz + (ix - 5)*dimy*dimz - 5];
                    // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
                    for (int iz = 5; iz < dimz + 5; iz++)
                    {
                        B[iz] = t0 * A[iz] +
                                t1z * (A[iz + 1] + A[iz - 1]) +
                                t2z * (A[iz + 2] + A[iz - 2]) +
                                t3z * (A[iz + 3] + A[iz - 3]) +
                                t4z * (A[iz + 4] + A[iz - 4]) +
                                t5z * (A[iz + 5] + A[iz - 5]);
                    }

                    for (int iz = 5; iz < dimz + 5; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + iys] + A[iz - iys]) +
                                t2x * (A[iz + 2*iys] + A[iz - 2*iys]) +
                                t3x * (A[iz + 3*iys] + A[iz - 3*iys]) +
                                t4x * (A[iz + 4*iys] + A[iz - 4*iys]) +
                                t5x * (A[iz + 5*iys] + A[iz - 5*iys]);
                    }

                    for (int iz = 5; iz < dimz + 5; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + ixs] + A[iz - ixs]) +
                                t2x * (A[iz + 2*ixs] + A[iz - 2*ixs]) +
                                t3x * (A[iz + 3*ixs] + A[iz - 3*ixs]) +
                                t4x * (A[iz + 4*ixs] + A[iz - 4*ixs]) +
                                t5x * (A[iz + 5*ixs] + A[iz - 5*ixs]);
                    }                   /* end for */

                    for (int iz = 5; iz < dimz + 5; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + ixs + iys] + A[iz - ixs - iys]) +
                                t2x * (A[iz + 2*ixs + 2*iys] + A[iz - 2*ixs - 2*iys]) +
                                t3x * (A[iz + 3*ixs + 3*iys] + A[iz - 3*ixs - 3*iys]) +
                                t4x * (A[iz + 4*ixs + 4*iys] + A[iz - 4*ixs - 4*iys]) +
                                t5x * (A[iz + 5*ixs + 5*iys] + A[iz - 5*ixs - 5*iys]);
                    }                   /* end for */

                }
            }

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }


    /* Return the diagonal component of the operator */
    return (double)std::real(t0);


}                               /* end app10_del2 */



    template <typename RmgType>
double FiniteDiff::app_cil_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    RmgTimer RT("App_cil: computation");
    int incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    RmgType ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    RmgType ihx, ihy, ihz;
    RmgType a1, a2, a3;
    RmgType ONE_t = 1.0;
    RmgType TWO_t = 2.0;
    RmgType THREE_t = 3.0;
    RmgType FOUR_t = 4.0;
    RmgType FIVE_t = 5.0;
    RmgType SIX_t = 6.0;
    RmgType EIGHT_t = 8.0;
    RmgType NINE_t = 9.0;
    RmgType TWELVE_t = 12.0;
    RmgType EIGHTTEEN_t = 18.0;
    RmgType TWENTYFOUR_t = 24.0;
    RmgType THIRTYFOUR_t = 34.0;
    RmgType THIRTYSIX_t = 36.0;
    RmgType FORTYEIGHT_t = 48.0;

    int ibrav = L->get_ibrav_type();


    ihx = 1.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    ihy = 1.0 / (gridhy * gridhy * L->get_yside() * L->get_yside());
    ihz = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    switch(ibrav) {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            if (FiniteDiff::check_anisotropy(gridhx, gridhy, gridhz, 0.0000001))
            {

                cc = (-FOUR_t / THREE_t) * (ihx + ihx + ihx);
                fcx = (FIVE_t / SIX_t) * ihx + (cc / EIGHT_t);
                ecxy = (ONE_t / TWELVE_t) * (ihx + ihx);
                incy = dimz + 2;
                incx = (dimz + 2) * (dimy + 2);
                incyr = dimz;
                incxr = dimz * dimy;

                for (int ix = 1; ix <= dimx; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;
                    for (int iy = 1; iy <= dimy; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (int iz = 1; iz <= dimz; iz++)
                        {

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                cc * rptr[ixs + iys + iz] +
                                fcx * (rptr[ixms + iys + iz] +
                                        rptr[ixps + iys + iz] +
                                        rptr[ixs + iyms + iz] +
                                        rptr[ixs + iyps + iz] +
                                        rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                ecxy * (rptr[ixms + iys + iz - 1] +
                                        rptr[ixps + iys + iz - 1] +
                                        rptr[ixs + iyms + iz - 1] +
                                        rptr[ixs + iyps + iz - 1] +
                                        rptr[ixms + iyms + iz] +
                                        rptr[ixms + iyps + iz] +
                                        rptr[ixps + iyms + iz] +
                                        rptr[ixps + iyps + iz] +
                                        rptr[ixms + iys + iz + 1] +
                                        rptr[ixps + iys + iz + 1] +
                                        rptr[ixs + iyms + iz + 1] + rptr[ixs + iyps + iz + 1]);


                        }           /* end for */

                    }               /* end for */

                }                   /* end for */
            }
            else
            {

                /* Compute coefficients for this grid spacing */
                cc = (-FOUR_t / THREE_t) * (ihx + ihy + ihz);

                fcx = FIVE_t/SIX_t * ihx + (cc / EIGHT_t);
                fcy = FIVE_t/SIX_t * ihy + (cc / EIGHT_t);
                fcz = FIVE_t/SIX_t * ihz + (cc / EIGHT_t);

                ecxy = (ONE_t / TWELVE_t) * (ihx + ihy);
                ecxz = (ONE_t / TWELVE_t) * (ihx + ihz);
                ecyz = (ONE_t / TWELVE_t) * (ihy + ihz);


                incy = dimz + 2;
                incx = (dimz + 2) * (dimy + 2);
                incyr = dimz;
                incxr = dimz * dimy;



                for (int ix = 1; ix <= dimx; ix++)
                {
                    ixs = ix * incx;
                    ixms = (ix - 1) * incx;
                    ixps = (ix + 1) * incx;
                    for (int iy = 1; iy <= dimy; iy++)
                    {
                        iys = iy * incy;
                        iyms = (iy - 1) * incy;
                        iyps = (iy + 1) * incy;

                        for (int iz = 1; iz <= dimz; iz++)
                        {

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                cc * rptr[ixs + iys + iz] +
                                fcx * rptr[ixms + iys + iz] +
                                fcx * rptr[ixps + iys + iz] +
                                fcy * rptr[ixs + iyms + iz] +
                                fcy * rptr[ixs + iyps + iz] +
                                fcz * rptr[ixs + iys + (iz - 1)] + fcz * rptr[ixs + iys + (iz + 1)];

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                ecxz * rptr[ixms + iys + iz - 1] +
                                ecxz * rptr[ixps + iys + iz - 1] +
                                ecyz * rptr[ixs + iyms + iz - 1] +
                                ecyz * rptr[ixs + iyps + iz - 1] +
                                ecxy * rptr[ixms + iyms + iz] +
                                ecxy * rptr[ixms + iyps + iz] +
                                ecxy * rptr[ixps + iyms + iz] +
                                ecxy * rptr[ixps + iyps + iz] +
                                ecxz * rptr[ixms + iys + iz + 1] +
                                ecxz * rptr[ixps + iys + iz + 1] +
                                ecyz * rptr[ixs + iyms + iz + 1] + ecyz * rptr[ixs + iyps + iz + 1];


                        }           /* end for */

                    }               /* end for */

                }                   /* end for */

            }                       /* end if */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            cc = ((-THREE_t / FOUR_t) * ihz) - ((FIVE_t / THREE_t) * ihx);
            a1 = ((THREE_t / EIGHT_t) * ihz) - ((ONE_t / SIX_t) * ihx);
            a2 = ((FIVE_t / EIGHTTEEN_t) * ihx) - ((ONE_t / TWENTYFOUR_t) * ihz);
            a3 = ((ONE_t / FORTYEIGHT_t) * ihz) + ((ONE_t / THIRTYSIX_t) * ihx);
            cc = TWO_t * cc;
            a1 = TWO_t * a1;
            a2 = TWO_t * a2;
            a3 = TWO_t * a3;

            for (int ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;
                for (int iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = 1; iz <= dimz; iz++)
                    {

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                            cc * rptr[ixs + iys + iz] +
                            a3 * (rptr[ixps + iys + iz - 1] +
                                    rptr[ixps + iyms + iz - 1] +
                                    rptr[ixs + iyms + iz - 1] +
                                    rptr[ixms + iys + iz - 1] +
                                    rptr[ixms + iyps + iz - 1] +
                                    rptr[ixs + iyps + iz - 1] +
                                    rptr[ixps + iys + iz + 1] +
                                    rptr[ixps + iyms + iz + 1] +
                                    rptr[ixs + iyms + iz + 1] +
                                    rptr[ixms + iys + iz + 1] +
                                    rptr[ixms + iyps + iz + 1] + 
                                    rptr[ixs + iyps + iz + 1]);


                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                            a2 * (rptr[ixps + iys + iz] +
                                    rptr[ixps + iyms + iz] +
                                    rptr[ixs + iyms + iz] +
                                    rptr[ixms + iys + iz] + 
                                    rptr[ixms + iyps + iz] + 
                                    rptr[ixs + iyps + iz]);

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                            a1 * (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 1]);

                    }               /* end for */

                }                   /* end for */

            }                       /* end for */

            break;

        case CUBIC_FC:

            cc = -(THIRTYFOUR_t / SIX_t) * ihx;
            a1 = (FOUR_t / NINE_t) * ihx;
            a2 = (ONE_t / EIGHTTEEN_t) * ihx;

            for (int ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;
                for (int iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = 1; iz <= dimz; iz++)
                    {

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                            cc * rptr[ix * incx + iys + iz] +
                            a1 * rptr[ixms + iys + iz] +
                            a1 * rptr[ixms + iys + iz + 1] +
                            a1 * rptr[ixms + iyps + iz] +
                            a1 * rptr[ixs + iyms + iz] +
                            a1 * rptr[ixs + iyms + iz + 1] +
                            a1 * rptr[ixs + iys + iz - 1] +
                            a1 * rptr[ixs + iys + iz + 1] +
                            a1 * rptr[ixs + iyps + iz - 1] +
                            a1 * rptr[ixs + iyps + iz] +
                            a1 * rptr[ixps + iyms + iz] +
                            a1 * rptr[ixps + iys + iz - 1] + 
                            a1 * rptr[ixps + iys + iz];


                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                            a2 * rptr[ixms + iyms + iz + 1] +
                            a2 * rptr[ixms + iyps + iz - 1] +
                            a2 * rptr[ixms + iyps + iz + 1] +
                            a2 * rptr[ixps + iyms + iz - 1] +
                            a2 * rptr[ixps + iyms + iz + 1] + 
                            a2 * rptr[ixps + iyps + iz - 1];

                    }               /* end for */

                }                   /* end for */

            }                       /* end for */
            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Grid symmetry not programmed yet in app_cil_fourth.\n");

    } // end switch

    return (double)std::real(cc);

}                               /* end app_cil */


    template <typename RmgType>
void FiniteDiff::app_cir_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz)
{

    RmgTimer RT("App_cir: computation");
    int ixs, iys, ixms, ixps, iyms, iyps;
    RmgType c000(0.5), c100(1.0/12.0);
    RmgType Bc(2.0 / 3.0);
    RmgType Bf(1.0 / 36.0);
    RmgType Bch(7.0 / 12.0);
    RmgType Bfh(1.0 / 24.0);
    RmgType Bz(1.0 / 12.0);

    int ibrav = L->get_ibrav_type();
    int incx = (dimz + 2) * (dimy + 2);
    int incy = dimz + 2;
    int incxr = dimz * dimy;
    int incyr = dimz;

    switch(ibrav) {

        case CUBIC_PRIMITIVE: 
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (int ix = 1; ix < dimx + 1; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;

                for (int iy = 1; iy < dimy + 1; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;

                    for (int iz = 1; iz < dimz + 1; iz++)
                    {

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                            c100 * (rptr[ixs + iys + (iz - 1)] +
                                    rptr[ixs + iys + (iz + 1)] +
                                    rptr[ixms + iys + iz] +
                                    rptr[ixps + iys + iz] +
                                    rptr[ixs + iyms + iz] +
                                    rptr[ixs + iyps + iz]) + 
                            c000 *  rptr[ixs + iys + iz];

                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            for(int ix = 1;ix < dimx + 1;ix++) {

                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;
                for(int iy = 1;iy < dimy + 1;iy++) {

                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;
                    for(int iz = 1;iz < dimz + 1;iz++) {

                        b[(ix - 1)*incxr + (iy - 1)*incyr + (iz - 1)] =
                            Bch * rptr[ixs + iys + iz] +
                            Bz * rptr[ixs + iys + iz - 1] +
                            Bz * rptr[ixs + iys + iz + 1] +
                            Bfh * rptr[ixps + iys + iz] +
                            Bfh * rptr[ixps + iyms + iz] +
                            Bfh * rptr[ixs + iyms + iz] +
                            Bfh * rptr[ixms + iys + iz] +
                            Bfh * rptr[ixms + iyps + iz] +
                            Bfh * rptr[ixs + iyps + iz];

                    } /* end for */

                } /* end for */

            } /* end for */
            break;

        case CUBIC_FC:

            for (int ix = 1; ix <= dimx; ix++)
            {
                ixs = ix * incx;
                ixms = (ix - 1) * incx;
                ixps = (ix + 1) * incx;
                for (int iy = 1; iy <= dimy; iy++)
                {
                    iys = iy * incy;
                    iyms = (iy - 1) * incy;
                    iyps = (iy + 1) * incy;
                    for (int iz = 1; iz <= dimz; iz++)
                    {

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                            Bc * rptr[ixs + iys + iz] +
                            Bf * (rptr[ixms + iys + iz] +
                                    rptr[ixms + iys + iz + 1] +
                                    rptr[ixms + iyps + iz] +
                                    rptr[ixs + iyms + iz] +
                                    rptr[ixs + iyms + iz + 1] +
                                    rptr[ixs + iys + iz - 1] +
                                    rptr[ixs + iys + iz + 1] +
                                    rptr[ixs + iyps + iz - 1] +
                                    rptr[ixs + iyps + iz] +
                                    rptr[ixps + iyms + iz] + 
                                    rptr[ixps + iys + iz - 1] + 
                                    rptr[ixps + iys + iz]);


                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Grid symmetry not programmed yet in app_cir_fourth.\n");

    }

}


    template <typename RmgType>
void FiniteDiff::app_gradient_sixth (RmgType * rptr, RmgType * wxr, RmgType *wyr, RmgType *wzr, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{

    int ixs = (dimy + 6) * (dimz + 6);
    int iys = (dimz + 6);
    int ix1 = dimy * dimz;
    int iy1 = dimz;

    int ibrav = L->get_ibrav_type();

    RmgType t1x (3.0 / ( 4.0 * gridhx * L->get_xside()));
    RmgType t2x (-3.0 / (20.0 * gridhx * L->get_xside()));
    RmgType t3x (1.0 / (60.0 * gridhx * L->get_xside()));

    RmgType t1y (3.0 / ( 4.0 * gridhy * L->get_yside()));
    RmgType t2y (-3.0 / (20.0 * gridhy * L->get_yside()));
    RmgType t3y (1.0 / (60.0 * gridhy * L->get_yside()));

    RmgType t1z (3.0 / ( 4.0 * gridhz * L->get_zside()));
    RmgType t2z (-3.0 / (20.0 * gridhz * L->get_zside()));
    RmgType t3z (1.0 / (60.0 * gridhz * L->get_zside()));

    switch (ibrav)
    {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (int ix = 3; ix < dimx + 3; ix++)
            {

                for (int iy = 3; iy < dimy + 3; iy++)
                {

                    for (int iz = 3; iz < dimz + 3; iz++)
                    {

                        wxr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                            -t3x * rptr[(ix - 3) * ixs + iy * iys + iz] +
                            -t2x * rptr[(ix - 2) * ixs + iy * iys + iz] +
                            -t1x * rptr[(ix - 1) * ixs + iy * iys + iz] +
                            t1x * rptr[(ix + 1) * ixs + iy * iys + iz] +
                            t2x * rptr[(ix + 2) * ixs + iy * iys + iz] +
                            t3x * rptr[(ix + 3) * ixs + iy * iys + iz];

                        wyr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                            -t3y * rptr[ix * ixs + (iy - 3) * iys + iz] +
                            -t2y * rptr[ix * ixs + (iy - 2) * iys + iz] +
                            -t1y * rptr[ix * ixs + (iy - 1) * iys + iz] +
                            t1y * rptr[ix * ixs + (iy + 1) * iys + iz] +
                            t2y * rptr[ix * ixs + (iy + 2) * iys + iz] +
                            t3y * rptr[ix * ixs + (iy + 3) * iys + iz];

                        wzr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                            -t3z * rptr[ix * ixs + iy * iys + iz - 3] +
                            -t2z * rptr[ix * ixs + iy * iys + iz - 2] +
                            -t1z * rptr[ix * ixs + iy * iys + iz - 1] +
                            t1z * rptr[ix * ixs + iy * iys + iz + 1] +
                            t2z * rptr[ix * ixs + iy * iys + iz + 2] +
                            t3z * rptr[ix * ixs + iy * iys + iz + 3];

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }                           /* end switch */
}

    template <typename RmgType>
void FiniteDiff::app_gradient_eighth (RmgType * __restrict__ rptr, RmgType * __restrict__ wxr, RmgType * __restrict__ wyr, RmgType * __restrict__ wzr, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{

    int ixs = (dimy + 8) * (dimz + 8);
    int iys = (dimz + 8);
    int ix1 = dimy * dimz;
    int iy1 = dimz;

    int ibrav = L->get_ibrav_type();

    RmgType t1x (4.0 / ( 5.0 * gridhx * L->get_xside()));
    RmgType t2x (-1.0 / (5.0 * gridhx * L->get_xside()));
    RmgType t3x (4.0 / (105.0 * gridhx * L->get_xside()));
    RmgType t4x (-1.0 / (280.0 * gridhx * L->get_xside()));

    RmgType t1y (4.0 / ( 5.0 * gridhy * L->get_yside()));
    RmgType t2y (-1.0 / (5.0 * gridhy * L->get_yside()));
    RmgType t3y (4.0 / (105.0 * gridhy * L->get_yside()));
    RmgType t4y (-1.0 / (280.0 * gridhy * L->get_yside()));

    RmgType t1z (4.0/ ( 5.0 * gridhz * L->get_zside()));
    RmgType t2z (-1.0 / (5.0 * gridhz * L->get_zside()));
    RmgType t3z (4.0 / (105.0 * gridhz * L->get_zside()));
    RmgType t4z (-1.0 / (280.0 * gridhz * L->get_zside()));
    RmgType hex_t(0.5*1.154700538379);

    RmgType x_gx[4] ,x_gy[4] ,x_gz[4] ,y_gx[4] ,y_gy[4] ,y_gz[4] ,xy_gx[4] ,xy_gy[4] ,xy_gz[4] ,nxy_gx[4] ,nxy_gy[4] ,nxy_gz[4];

    for(int i = 0; i < 4; i++)
    {
        x_gx[i]   = (RmgType)LC->axis_x_gx[i]; 
        x_gy[i]   = (RmgType)LC->axis_x_gy[i]; 
        x_gz[i]   = (RmgType)LC->axis_x_gz[i]; 
        y_gx[i]   = (RmgType)LC->axis_y_gx[i]; 
        y_gy[i]   = (RmgType)LC->axis_y_gy[i]; 
        y_gz[i]   = (RmgType)LC->axis_y_gz[i]; 
        xy_gx[i]  = (RmgType)LC->axis_xy_gx[i];
        xy_gy[i]  = (RmgType)LC->axis_xy_gy[i];
        xy_gz[i]  = (RmgType)LC->axis_xy_gz[i];
        nxy_gx[i] = (RmgType)LC->axis_nxy_gx[i];
        nxy_gy[i] = (RmgType)LC->axis_nxy_gy[i];
        nxy_gz[i] = (RmgType)LC->axis_nxy_gz[i];
    }

    int id = 1;
    switch (ibrav)
    {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (int ix = 4; ix < dimx + 4; ix++)
            {

                for (int iy = 4; iy < dimy + 4; iy++)
                {

                    RmgType *A = &wxr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    RmgType *B = &rptr[ix * ixs + iy * iys];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4x * ( B[iz + 4*ixs] - B[iz - 4*ixs]) +
                            t3x * ( B[iz + 3*ixs] - B[iz - 3*ixs]) +
                            t2x * ( B[iz + 2*ixs] - B[iz - 2*ixs]) +
                            t1x * ( B[iz + ixs] - B[iz - ixs]);
                    }
                }                   /* end for */
            }                       /* end for */

            for (int ix = 4; ix < dimx + 4; ix++)
            {

                for (int iy = 4; iy < dimy + 4; iy++)
                {

                    RmgType *A = &wyr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    RmgType *B = &rptr[ix * ixs + iy * iys];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4y * ( B[iz + 4*iys] - B[iz - 4*iys]) +
                            t3y * ( B[iz + 3*iys] - B[iz - 3*iys]) +
                            t2y * ( B[iz + 2*iys] - B[iz - 2*iys]) +
                            t1y * ( B[iz + iys] - B[iz - iys]);

                    }
                }                   /* end for */
            }                       /* end for */

            for (int ix = 4; ix < dimx + 4; ix++)
            {

                for (int iy = 4; iy < dimy + 4; iy++)
                {

                    RmgType *A = &wzr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    RmgType *B = &rptr[ix * ixs + iy * iys];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4z * ( B[iz + 4] - B[iz - 4]) +
                            t3z * ( B[iz + 3] - B[iz - 3]) +
                            t2z * ( B[iz + 2] - B[iz - 2]) +
                            t1z * ( B[iz + 1] - B[iz - 1]);

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

            break;

        case HEXAGONAL2:
            id = -1;
        case HEXAGONAL:

            for (int ix = 4; ix < dimx + 4; ix++)
            {

                for (int iy = 4; iy < dimy + 4; iy++)
                {

                    RmgType *A = &wxr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    RmgType *B = &rptr[ix * ixs + iy * iys];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4x * ( B[iz + 4*ixs] - B[iz - 4*ixs]) +
                            t3x * ( B[iz + 3*ixs] - B[iz - 3*ixs]) +
                            t2x * ( B[iz + 2*ixs] - B[iz - 2*ixs]) +
                            t1x * ( B[iz + ixs] - B[iz - ixs]);
                    }

                    A = &wyr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            hex_t * t4y * ( B[iz + 4*iys] - B[iz - 4*iys]) +
                            hex_t * t3y * ( B[iz + 3*iys] - B[iz - 3*iys]) +
                            hex_t * t2y * ( B[iz + 2*iys] - B[iz - 2*iys]) +
                            hex_t * t1y * ( B[iz + iys] - B[iz - iys]) +

                            hex_t * t4y * ( -B[iz - id*4*ixs - 4*iys] + B[iz + id*4*ixs + 4*iys]) +
                            hex_t * t3y * ( -B[iz - id*3*ixs - 3*iys] + B[iz + id*3*ixs + 3*iys]) +
                            hex_t * t2y * ( -B[iz - id*2*ixs - 2*iys] + B[iz + id*2*ixs + 2*iys]) +
                            hex_t * t1y * ( -B[iz - id*ixs - iys] + B[iz + id*ixs + iys]);
                    }

                    A = &wzr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4z * ( B[iz + 4] - B[iz - 4]) +
                            t3z * ( B[iz + 3] - B[iz - 3]) +
                            t2z * ( B[iz + 2] - B[iz - 2]) +
                            t1z * ( B[iz + 1] - B[iz - 1]);

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

            break;

        case MONOCLINIC_PRIMITIVE:

            for (int ix = 4; ix < dimx + 4; ix++)
            {

                for (int iy = 4; iy < dimy + 4; iy++)
                {

                    RmgType *A = &wxr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    RmgType *B = &rptr[ix * ixs + iy * iys];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            x_gx[3] * ( B[iz + 4*ixs] - B[iz - 4*ixs]) +
                            x_gx[2] * ( B[iz + 3*ixs] - B[iz - 3*ixs]) +
                            x_gx[1] * ( B[iz + 2*ixs] - B[iz - 2*ixs]) +
                            x_gx[0] * ( B[iz + ixs] - B[iz - ixs]) +
                            y_gx[3] * ( B[iz + 4*iys] - B[iz - 4*iys]) +
                            y_gx[2] * ( B[iz + 3*iys] - B[iz - 3*iys]) +
                            y_gx[1] * ( B[iz + 2*iys] - B[iz - 2*iys]) +
                            y_gx[0] * ( B[iz + iys] - B[iz - iys]) +
                            xy_gx[3] * ( -B[iz - 4*ixs - 4*iys] + B[iz + 4*ixs + 4*iys]) +
                            xy_gx[2] * ( -B[iz - 3*ixs - 3*iys] + B[iz + 3*ixs + 3*iys]) +
                            xy_gx[1] * ( -B[iz - 2*ixs - 2*iys] + B[iz + 2*ixs + 2*iys]) +
                            xy_gx[0] * ( -B[iz - ixs - iys] + B[iz + ixs + iys]) +
                            nxy_gx[3] * ( -B[iz + 4*ixs - 4*iys] + B[iz - 4*ixs + 4*iys]) +
                            nxy_gx[2] * ( -B[iz + 3*ixs - 3*iys] + B[iz - 3*ixs + 3*iys]) +
                            nxy_gx[1] * ( -B[iz + 2*ixs - 2*iys] + B[iz - 2*ixs + 2*iys]) +
                            nxy_gx[0] * ( -B[iz + ixs - iys] + B[iz - ixs + iys]);


                    }

                    A = &wyr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            x_gy[3] * ( B[iz + 4*ixs] - B[iz - 4*ixs]) +
                            x_gy[2] * ( B[iz + 3*ixs] - B[iz - 3*ixs]) +
                            x_gy[1] * ( B[iz + 2*ixs] - B[iz - 2*ixs]) +
                            x_gy[0] * ( B[iz + ixs] - B[iz - ixs]) +
                            y_gy[3] * ( B[iz + 4*iys] - B[iz - 4*iys]) +
                            y_gy[2] * ( B[iz + 3*iys] - B[iz - 3*iys]) +
                            y_gy[1] * ( B[iz + 2*iys] - B[iz - 2*iys]) +
                            y_gy[0] * ( B[iz + iys] - B[iz - iys]) +
                            xy_gy[3] * ( -B[iz - 4*ixs - 4*iys] + B[iz + 4*ixs + 4*iys]) +
                            xy_gy[2] * ( -B[iz - 3*ixs - 3*iys] + B[iz + 3*ixs + 3*iys]) +
                            xy_gy[1] * ( -B[iz - 2*ixs - 2*iys] + B[iz + 2*ixs + 2*iys]) +
                            xy_gy[0] * ( -B[iz - ixs - iys] + B[iz + ixs + iys]) +
                            nxy_gy[3] * ( -B[iz + 4*ixs - 4*iys] + B[iz - 4*ixs + 4*iys]) +
                            nxy_gy[2] * ( -B[iz + 3*ixs - 3*iys] + B[iz - 3*ixs + 3*iys]) +
                            nxy_gy[1] * ( -B[iz + 2*ixs - 2*iys] + B[iz - 2*ixs + 2*iys]) +
                            nxy_gy[0] * ( -B[iz + ixs - iys] + B[iz - ixs + iys]);
                    }

                    A = &wzr[(ix - 4) * ix1 + (iy - 4) * iy1 - 4];
                    for (int iz = 4; iz < dimz + 4; iz++)
                    {
                        A[iz] =
                            t4z * ( B[iz + 4] - B[iz - 4]) +
                            t3z * ( B[iz + 3] - B[iz - 3]) +
                            t2z * ( B[iz + 2] - B[iz - 2]) +
                            t1z * ( B[iz + 1] - B[iz - 1]);
                    }                   /* end for */
                }
            }

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }                           /* end switch */


}


    template <typename RmgType>
void FiniteDiff::app_gradient_tenth (RmgType * __restrict__ rptr, RmgType * __restrict__ wxr, RmgType * __restrict__ wyr, RmgType * __restrict__ wzr, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{

    int ixs = (dimy + 10) * (dimz + 10);
    int iys = (dimz + 10);
    int ix1 = dimy * dimz;
    int iy1 = dimz;

    int ibrav = L->get_ibrav_type();

    RmgType t1x (0.833333333333333370 / (gridhx * L->get_xside()));
    RmgType t2x (-0.238095238095238082 / (gridhx * L->get_xside()));
    RmgType t3x (0.059523809523809521 / (gridhx * L->get_xside()));
    RmgType t4x (-0.009920634920634920 / (gridhx * L->get_xside()));
    RmgType t5x (0.000793650793650794 / (gridhx * L->get_xside()));

    RmgType t1y (0.833333333333333370 / (gridhy * L->get_yside()));
    RmgType t2y (-0.238095238095238082 / (gridhy * L->get_yside()));
    RmgType t3y (0.059523809523809521 / (gridhy * L->get_yside()));
    RmgType t4y (-0.009920634920634920 / (gridhy * L->get_yside()));
    RmgType t5y (0.000793650793650794 / (gridhy * L->get_yside()));

    RmgType t1z (0.833333333333333370/ (gridhz * L->get_zside()));
    RmgType t2z (-0.238095238095238082 / (gridhz * L->get_zside()));
    RmgType t3z (0.059523809523809521 / (gridhz * L->get_zside()));
    RmgType t4z (-0.009920634920634920 / (gridhz * L->get_zside()));
    RmgType t5z (0.000793650793650794 / (gridhz * L->get_zside()));


    switch (ibrav)
    {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (int ix = 5; ix < dimx + 5; ix++)
            {

                for (int iy = 5; iy < dimy + 5; iy++)
                {

                    for (int iz = 5; iz < dimz + 5; iz++)
                    {

                        wxr[(ix - 5) * ix1 + (iy - 5) * iy1 + iz - 5] =
                            -t5x * rptr[(ix - 5) * ixs + iy * iys + iz] +
                            -t4x * rptr[(ix - 4) * ixs + iy * iys + iz] +
                            -t3x * rptr[(ix - 3) * ixs + iy * iys + iz] +
                            -t2x * rptr[(ix - 2) * ixs + iy * iys + iz] +
                            -t1x * rptr[(ix - 1) * ixs + iy * iys + iz] +
                            t1x * rptr[(ix + 1) * ixs + iy * iys + iz] +
                            t2x * rptr[(ix + 2) * ixs + iy * iys + iz] +
                            t3x * rptr[(ix + 3) * ixs + iy * iys + iz] +
                            t4x * rptr[(ix + 4) * ixs + iy * iys + iz] +
                            t5x * rptr[(ix + 5) * ixs + iy * iys + iz];

                        wyr[(ix - 5) * ix1 + (iy - 5) * iy1 + iz - 5] =
                            -t5y * rptr[ix * ixs + (iy - 5) * iys + iz] +
                            -t4y * rptr[ix * ixs + (iy - 4) * iys + iz] +
                            -t3y * rptr[ix * ixs + (iy - 3) * iys + iz] +
                            -t2y * rptr[ix * ixs + (iy - 2) * iys + iz] +
                            -t1y * rptr[ix * ixs + (iy - 1) * iys + iz] +
                            t1y * rptr[ix * ixs + (iy + 1) * iys + iz] +
                            t2y * rptr[ix * ixs + (iy + 2) * iys + iz] +
                            t3y * rptr[ix * ixs + (iy + 3) * iys + iz] +
                            t4y * rptr[ix * ixs + (iy + 4) * iys + iz] +
                            t5y * rptr[ix * ixs + (iy + 5) * iys + iz];

                        wzr[(ix - 5) * ix1 + (iy - 5) * iy1 + iz - 5] =
                            -t5z * rptr[ix * ixs + iy * iys + iz - 5] +
                            -t4z * rptr[ix * ixs + iy * iys + iz - 4] +
                            -t3z * rptr[ix * ixs + iy * iys + iz - 3] +
                            -t2z * rptr[ix * ixs + iy * iys + iz - 2] +
                            -t1z * rptr[ix * ixs + iy * iys + iz - 1] +
                            t1z * rptr[ix * ixs + iy * iys + iz + 1] +
                            t2z * rptr[ix * ixs + iy * iys + iz + 2] +
                            t3z * rptr[ix * ixs + iy * iys + iz + 3] +
                            t4z * rptr[ix * ixs + iy * iys + iz + 4] +
                            t5z * rptr[ix * ixs + iy * iys + iz + 5];


                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }                           /* end switch */


}



// Openmp version. Very simple with no cache optimizations as of yet.
    template <typename RmgType>
double FiniteDiff::app_cil_fourth_threaded (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz)
{

    RmgTimer RT("App_cil: computation");
    int incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    RmgType ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    RmgType ihx, ihy, ihz;
    RmgType a1, a2, a3;
    RmgType ONE_t = 1.0;
    RmgType TWO_t = 2.0;
    RmgType THREE_t = 3.0;
    RmgType FOUR_t = 4.0;
    RmgType FIVE_t = 5.0;
    RmgType SIX_t = 6.0;
    RmgType EIGHT_t = 8.0;
    RmgType NINE_t = 9.0;
    RmgType TWELVE_t = 12.0;
    RmgType EIGHTTEEN_t = 18.0;
    RmgType TWENTYFOUR_t = 24.0;
    RmgType THIRTYFOUR_t = 34.0;
    RmgType THIRTYSIX_t = 36.0;
    RmgType FORTYEIGHT_t = 48.0;

    int ibrav = L->get_ibrav_type();


    ihx = 1.0 / (gridhx * gridhx * L->get_xside() * L->get_xside());
    ihy = 1.0 / (gridhy * gridhy * L->get_yside() * L->get_yside());
    ihz = 1.0 / (gridhz * gridhz * L->get_zside() * L->get_zside());


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    switch(ibrav) {

        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            /* Compute coefficients for this grid spacing */
            cc = (-FOUR_t / THREE_t) * (ihx + ihy + ihz);

            fcx = FIVE_t/SIX_t * ihx + (cc / EIGHT_t);
            fcy = FIVE_t/SIX_t * ihy + (cc / EIGHT_t);
            fcz = FIVE_t/SIX_t * ihz + (cc / EIGHT_t);

            ecxy = (ONE_t / TWELVE_t) * (ihx + ihy);
            ecxz = (ONE_t / TWELVE_t) * (ihx + ihz);
            ecyz = (ONE_t / TWELVE_t) * (ihy + ihz);


            incy = dimz + 2;
            incx = (dimz + 2) * (dimy + 2);
            incyr = dimz;
            incxr = dimz * dimy;

            int ix, iy, iz;
#pragma omp parallel private(ix,iy,iz,ixs,ixms,ixps,iys,iyms,iyps)
            {
#pragma omp for schedule(static, 2) nowait
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

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                cc * rptr[ixs + iys + iz] +
                                fcx * rptr[ixms + iys + iz] +
                                fcx * rptr[ixps + iys + iz] +
                                fcy * rptr[ixs + iyms + iz] +
                                fcy * rptr[ixs + iyps + iz] +
                                fcz * rptr[ixs + iys + (iz - 1)] + fcz * rptr[ixs + iys + (iz + 1)];

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                ecxz * rptr[ixms + iys + iz - 1] +
                                ecxz * rptr[ixps + iys + iz - 1] +
                                ecyz * rptr[ixs + iyms + iz - 1] +
                                ecyz * rptr[ixs + iyps + iz - 1] +
                                ecxy * rptr[ixms + iyms + iz] +
                                ecxy * rptr[ixms + iyps + iz] +
                                ecxy * rptr[ixps + iyms + iz] +
                                ecxy * rptr[ixps + iyps + iz] +
                                ecxz * rptr[ixms + iys + iz + 1] +
                                ecxz * rptr[ixps + iys + iz + 1] +
                                ecyz * rptr[ixs + iyms + iz + 1] + ecyz * rptr[ixs + iyps + iz + 1];


                        }           /* end for */

                    }               /* end for */

                }                   /* end for */
            }
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            cc = ((-THREE_t / FOUR_t) * ihz) - ((FIVE_t / THREE_t) * ihx);
            a1 = ((THREE_t / EIGHT_t) * ihz) - ((ONE_t / SIX_t) * ihx);
            a2 = ((FIVE_t / EIGHTTEEN_t) * ihx) - ((ONE_t / TWENTYFOUR_t) * ihz);
            a3 = ((ONE_t / FORTYEIGHT_t) * ihz) + ((ONE_t / THIRTYSIX_t) * ihx);
            cc = TWO_t * cc;
            a1 = TWO_t * a1;
            a2 = TWO_t * a2;
            a3 = TWO_t * a3;

#pragma omp parallel private(ix,iy,iz,ixs,ixms,ixps,iys,iyms,iyps)
            {
#pragma omp for schedule(static, 2) nowait

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

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                cc * rptr[ixs + iys + iz] +
                                a3 * (rptr[ixps + iys + iz - 1] +
                                        rptr[ixps + iyms + iz - 1] +
                                        rptr[ixs + iyms + iz - 1] +
                                        rptr[ixms + iys + iz - 1] +
                                        rptr[ixms + iyps + iz - 1] +
                                        rptr[ixs + iyps + iz - 1] +
                                        rptr[ixps + iys + iz + 1] +
                                        rptr[ixps + iyms + iz + 1] +
                                        rptr[ixs + iyms + iz + 1] +
                                        rptr[ixms + iys + iz + 1] +
                                        rptr[ixms + iyps + iz + 1] + 
                                        rptr[ixs + iyps + iz + 1]);


                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                a2 * (rptr[ixps + iys + iz] +
                                        rptr[ixps + iyms + iz] +
                                        rptr[ixs + iyms + iz] +
                                        rptr[ixms + iys + iz] + 
                                        rptr[ixms + iyps + iz] + 
                                        rptr[ixs + iyps + iz]);

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                a1 * (rptr[ixs + iys + iz - 1] + rptr[ixs + iys + iz + 1]);

                        }               /* end for */

                    }                   /* end for */

                }                       /* end for */
            }
            break;

        case CUBIC_FC:

            cc = -(THIRTYFOUR_t / SIX_t) * ihx;
            a1 = (FOUR_t / NINE_t) * ihx;
            a2 = (ONE_t / EIGHTTEEN_t) * ihx;

#pragma omp parallel private(ix,iy,iz,ixs,ixms,ixps,iys,iyms,iyps)
            {
#pragma omp for schedule(static, 2) nowait

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

                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                                cc * rptr[ix * incx + iys + iz] +
                                a1 * rptr[ixms + iys + iz] +
                                a1 * rptr[ixms + iys + iz + 1] +
                                a1 * rptr[ixms + iyps + iz] +
                                a1 * rptr[ixs + iyms + iz] +
                                a1 * rptr[ixs + iyms + iz + 1] +
                                a1 * rptr[ixs + iys + iz - 1] +
                                a1 * rptr[ixs + iys + iz + 1] +
                                a1 * rptr[ixs + iyps + iz - 1] +
                                a1 * rptr[ixs + iyps + iz] +
                                a1 * rptr[ixps + iyms + iz] +
                                a1 * rptr[ixps + iys + iz - 1] + 
                                a1 * rptr[ixps + iys + iz];


                            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                                a2 * rptr[ixms + iyms + iz + 1] +
                                a2 * rptr[ixms + iyps + iz - 1] +
                                a2 * rptr[ixms + iyps + iz + 1] +
                                a2 * rptr[ixps + iyms + iz - 1] +
                                a2 * rptr[ixps + iyms + iz + 1] + 
                                a2 * rptr[ixps + iyps + iz - 1];

                        }               /* end for */

                    }                   /* end for */

                }                       /* end for */
            }
            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Grid symmetry not programmed yet in app_cil_fourth.\n");

    } // end switch

    return (double)std::real(cc);

}                               /* end app_cil */

#include "rmg_complex.h"


template <typename RmgType>
double FiniteDiff::app8_combined(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz, double *kvec)
{
    return FiniteDiff::app8_combined(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec, false);
}

template <typename RmgType>
double FiniteDiff::app8_combined(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz, double *kvec, bool use_gpu)
{
    int ibrav = L->get_ibrav_type();

//    if(ibrav == CUBIC_PRIMITIVE || ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
//    {
//        return FiniteDiff::app8_combined_orthorhombic(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec, use_gpu);
//    }

    return FiniteDiff::app8_combined_general(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, kvec, use_gpu);

} /* end app8_combined */



template <typename RmgType>
double FiniteDiff::app8_combined_orthorhombic(
		RmgType * __restrict__ a, RmgType * __restrict__ b, 
		int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
		double *kvec, bool use_gpu)
{
    double xside = L->get_xside();
    double yside = L->get_yside();
    double zside = L->get_zside();
    double s1 = 2.0;
    fdparms_o8<RmgType> c;

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

    int ixs = (dimy + 8) * (dimz + 8);
    int iys = (dimz + 8);

    // nine and seven point stencils, 2nd derivative, extrapolated
    int ic = 4;
    double x[10], w1[10], w2[10];
    for(int i=0;i<10;i++) x[i] = (double)i;
    gen_weights(9, 2, (double)ic, x, w1);
    gen_weights(7, 2, (double)(ic-1), x, w2);
    double hf = 1.0, c1, c2=0.0;
    double h2x = gridhx * gridhx * xside * xside;
    double h2y = gridhy * gridhy * yside * yside;
    double h2z = gridhz * gridhz * zside * zside;

    double d1 = 6.0/560.0;
    double d2 = -8.0/3150.0;
    double dr = d1 / d2;
    double k2 = PI*PI/8.0;


    double maxh = std::max(h2x, h2y);
    maxh = std::max(maxh, h2z);

    double hadj = sqrt(h2x / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    double th2 = (c1*w1[ic] - c2*w2[ic-1]) / h2x;

    RmgType t1x ((c1*w1[ic+1] - c2*w2[ic]) * hf / h2x);
    RmgType t2x ((c1*w1[ic+2] - c2*w2[ic+1]) * hf / h2x);
    RmgType t3x ((c1*w1[ic+3] - c2*w2[ic+2]) * hf / h2x);
    RmgType t4x (c1*w1[ic+4] * hf / h2x);

    hadj = sqrt(h2y / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) / h2y;

    RmgType t1y ((c1*w1[ic+1] - c2*w2[ic]) * hf / h2y);
    RmgType t2y ((c1*w1[ic+2] - c2*w2[ic+1]) * hf / h2y);
    RmgType t3y ((c1*w1[ic+3] - c2*w2[ic+2]) * hf / h2y);
    RmgType t4y (c1*w1[ic+4] * hf / h2y);

    hadj = sqrt(h2z / maxh);
    //if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hf*hadj/k2);
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) /  h2z;
    RmgType t1z ((c1*w1[ic+1] - c2*w2[ic]) / h2z);
    RmgType t2z ((c1*w1[ic+2] - c2*w2[ic+1]) / h2z);
    RmgType t3z ((c1*w1[ic+3] - c2*w2[ic+2]) / h2z);
    RmgType t4z (c1*w1[ic+4] / h2z);
    RmgType t0 (th2);

    // When kvec[i] != 0 this includes the gradient component
    c.a0 = t0;
    c.gpt1x = t1x + s1*kvec[0] * I_t * (4.0 / ( 5.0 * gridhx * xside));
    c.gpt2x = t2x + s1*kvec[0] * I_t * (-1.0 / (5.0 * gridhx * xside));
    c.gpt3x = t3x + s1*kvec[0] * I_t * (4.0 / (105.0 * gridhx * xside));
    c.gpt4x = t4x + s1*kvec[0] * I_t * (-1.0 / (280.0 * gridhx * xside));

    c.gmt1x = t1x - s1*kvec[0] * I_t * (4.0 / ( 5.0 * gridhx * xside));
    c.gmt2x = t2x - s1*kvec[0] * I_t * (-1.0 / (5.0 * gridhx * xside));
    c.gmt3x = t3x - s1*kvec[0] * I_t * (4.0 / (105.0 * gridhx * xside));
    c.gmt4x = t4x - s1*kvec[0] * I_t * (-1.0 / (280.0 * gridhx * xside));

    c.gpt1y = t1y + s1*kvec[1] * I_t * (4.0 / ( 5.0 * gridhy * yside));
    c.gpt2y = t2y + s1*kvec[1] * I_t * (-1.0 / (5.0 * gridhy * yside));
    c.gpt3y = t3y + s1*kvec[1] * I_t * (4.0 / (105.0 * gridhy * yside));
    c.gpt4y = t4y + s1*kvec[1] * I_t *  (-1.0 / (280.0 * gridhy * yside));

    c.gmt1y = t1y - s1*kvec[1] * I_t * (4.0 / ( 5.0 * gridhy * yside));
    c.gmt2y = t2y - s1*kvec[1] * I_t * (-1.0 / (5.0 * gridhy * yside));
    c.gmt3y = t3y - s1*kvec[1] * I_t * (4.0 / (105.0 * gridhy * yside));
    c.gmt4y = t4y - s1*kvec[1] * I_t *  (-1.0 / (280.0 * gridhy * yside));

    c.gpt1z = t1z + s1*kvec[2] * I_t * (4.0/ ( 5.0 * gridhz * zside));
    c.gpt2z = t2z + s1*kvec[2] * I_t * (-1.0 / (5.0 * gridhz * zside));
    c.gpt3z = t3z + s1*kvec[2] * I_t * (4.0 / (105.0 * gridhz * zside));
    c.gpt4z = t4z + s1*kvec[2] * I_t * (-1.0 / (280.0 * gridhz * zside));

    c.gmt1z = t1z - s1*kvec[2] * I_t * (4.0/ ( 5.0 * gridhz * zside));
    c.gmt2z = t2z - s1*kvec[2] * I_t * (-1.0 / (5.0 * gridhz * zside));
    c.gmt3z = t3z - s1*kvec[2] * I_t * (4.0 / (105.0 * gridhz * zside));
    c.gmt4z = t4z - s1*kvec[2] * I_t * (-1.0 / (280.0 * gridhz * zside));

    // NULL b means we just want the diagonal component.
    if(b == NULL) return (double)std::real(t0);

#if HIP_ENABLED || CUDA_ENABLED
    if(use_gpu)
    {
        /* Return the diagonal component of the operator */
        app8_del2_gpu(a, b, dimx, dimy, dimz, c);
        return (double)std::real(t0);
    }
#endif

    for (int ix = 4; ix < dimx + 4; ix++)
    {

        for (int iy = 4; iy < dimy + 4; iy++)
        {

            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];

            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] = t0 * A[iz] +
                    c.gpt1z * A[iz + 1] + c.gmt1z * A[iz - 1] +
                    c.gpt2z * A[iz + 2] + c.gmt2z * A[iz - 2] +
                    c.gpt3z * A[iz + 3] + c.gmt3z * A[iz - 3] +
                    c.gpt4z * A[iz + 4] + c.gmt4z * A[iz - 4];
            }

            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    c.gpt1y * A[iz + iys] + c.gmt1y * A[iz - iys] +
                    c.gpt2y * A[iz + 2*iys] + c.gmt2y * A[iz - 2*iys] +
                    c.gpt3y * A[iz + 3*iys] + c.gmt3y * A[iz - 3*iys] +
                    c.gpt4y * A[iz + 4*iys] + c.gmt4y * A[iz - 4*iys];
            }

            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    c.gpt1x * A[iz + ixs] + c.gmt1x * A[iz - ixs] +
                    c.gpt2x * A[iz + 2*ixs] + c.gmt2x * A[iz - 2*ixs] +
                    c.gpt3x * A[iz + 3*ixs] + c.gmt3x * A[iz - 3*ixs] +
                    c.gpt4x * A[iz + 4*ixs] + c.gmt4x * A[iz - 4*ixs];
            }                   /* end for */

        }                       /* end for */
    }                           /* end for */

    /* Return the diagonal component of the operator */
    return (double)std::real(t0);

} /* end app8_combined_orthorhombic */


template <typename RmgType>
double FiniteDiff::app8_combined_general(RmgType * __restrict__ a, RmgType * __restrict__ b, 
		int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
		double *kvec, bool use_gpu)
{
    double s1 = 2.0;
    fdparms_o8<RmgType> c;

    int ibrav = L->get_ibrav_type();

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

    int ixs = (dimy + 8) * (dimz + 8);
    int iys = (dimz + 8);

    // nine and seven point stencils, 2nd derivative, extrapolated
    double c1, c2;
    double d1 = 6.0/560.0;
    double d2 = -8.0/3150.0;
    double dr = d1 / d2;
    double k2 = PI*PI/8.0;
    double th2=0.0, maxh, hadj;

    RmgType t1xy, t2xy, t3xy, t4xy, t1xz, t2xz, t3xz, t4xz, t1yz, t2yz, t3yz, t4yz;
    RmgType t1nxy, t2nxy, t3nxy, t4nxy, t1nxz, t2nxz, t3nxz, t4nxz, t1nyz, t2nyz, t3nyz, t4nyz;
    RmgType t1x, t2x, t3x, t4x;
    RmgType t1y, t2y, t3y, t4y;
    RmgType t1z, t2z, t3z, t4z;


    // Need to setup 9 axes
    maxh = std::max(LC->plane_dist_x, LC->plane_dist_y);
    maxh = std::max(maxh, LC->plane_dist_z);
    maxh = std::max(maxh, LC->plane_dist_xy);
    maxh = std::max(maxh, LC->plane_dist_nxy);
    maxh = std::max(maxh, LC->plane_dist_xz);
    maxh = std::max(maxh, LC->plane_dist_nxz);
    maxh = std::max(maxh, LC->plane_dist_yz);
    maxh = std::max(maxh, LC->plane_dist_nyz);

    hadj = sqrt(LC->plane_dists[0] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[0] - c2*LC_6->plane_centers[0];
    t1x = c1*LC->axis_lc[0][3] - c2*LC_6->axis_lc[0][2];
    t2x = c1*LC->axis_lc[0][2] - c2*LC_6->axis_lc[0][1];
    t3x = c1*LC->axis_lc[0][1] - c2*LC_6->axis_lc[0][0];
    t4x = c1*LC->axis_lc[0][0];
    hadj = sqrt(LC->plane_dists[1] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[1] - c2*LC_6->plane_centers[1];
    t1y = c1*LC->axis_lc[1][3] - c2*LC_6->axis_lc[1][2];
    t2y = c1*LC->axis_lc[1][2] - c2*LC_6->axis_lc[1][1];
    t3y = c1*LC->axis_lc[1][1] - c2*LC_6->axis_lc[1][0];
    t4y = c1*LC->axis_lc[1][0];

    hadj = sqrt(LC->plane_dists[2] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[2] - c2*LC_6->plane_centers[2];
    t1z = c1*LC->axis_lc[2][3] - c2*LC_6->axis_lc[2][2];
    t2z = c1*LC->axis_lc[2][2] - c2*LC_6->axis_lc[2][1];
    t3z = c1*LC->axis_lc[2][1] - c2*LC_6->axis_lc[2][0];
    t4z = c1*LC->axis_lc[2][0];

    hadj = sqrt(LC->plane_dists[3] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[3] - c2*LC_6->plane_centers[3];
    t1xy = c1*LC->axis_lc[3][3] - c2*LC_6->axis_lc[3][2];
    t2xy = c1*LC->axis_lc[3][2] - c2*LC_6->axis_lc[3][1];
    t3xy = c1*LC->axis_lc[3][1] - c2*LC_6->axis_lc[3][0];
    t4xy = c1*LC->axis_lc[3][0];

    hadj = sqrt(LC->plane_dists[6] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[6] - c2*LC_6->plane_centers[6];
    t1nxy = c1*LC->axis_lc[6][3] - c2*LC_6->axis_lc[6][2];
    t2nxy = c1*LC->axis_lc[6][2] - c2*LC_6->axis_lc[6][1];
    t3nxy = c1*LC->axis_lc[6][1] - c2*LC_6->axis_lc[6][0];
    t4nxy = c1*LC->axis_lc[6][0];

    hadj = sqrt(LC->plane_dists[4] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[4] - c2*LC_6->plane_centers[4];
    t1xz = c1*LC->axis_lc[4][3] - c2*LC_6->axis_lc[4][2];
    t2xz = c1*LC->axis_lc[4][2] - c2*LC_6->axis_lc[4][1];
    t3xz = c1*LC->axis_lc[4][1] - c2*LC_6->axis_lc[4][0];
    t4xz = c1*LC->axis_lc[4][0];

    hadj = sqrt(LC->plane_dists[5] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[5] - c2*LC_6->plane_centers[5];
    t1yz = c1*LC->axis_lc[5][3] - c2*LC_6->axis_lc[5][2];
    t2yz = c1*LC->axis_lc[5][2] - c2*LC_6->axis_lc[5][1];
    t3yz = c1*LC->axis_lc[5][1] - c2*LC_6->axis_lc[5][0];
    t4yz = c1*LC->axis_lc[5][0];

    hadj = sqrt(LC->plane_dists[7] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[7] - c2*LC_6->plane_centers[7];
    t1nxz = c1*LC->axis_lc[7][3] - c2*LC_6->axis_lc[7][2];
    t2nxz = c1*LC->axis_lc[7][2] - c2*LC_6->axis_lc[7][1];
    t3nxz = c1*LC->axis_lc[7][1] - c2*LC_6->axis_lc[7][0];
    t4nxz = c1*LC->axis_lc[7][0];

    hadj = sqrt(LC->plane_dists[8] / maxh);
    c2 = 0.0;
    if(this->alt_laplacian) c2 = -1.0/(1.0+dr*hadj/k2);
if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += c1*LC->plane_centers[8] - c2*LC_6->plane_centers[8];
    t1nyz = c1*LC->axis_lc[8][3] - c2*LC_6->axis_lc[8][2];
    t2nyz = c1*LC->axis_lc[8][2] - c2*LC_6->axis_lc[8][1];
    t3nyz = c1*LC->axis_lc[8][1] - c2*LC_6->axis_lc[8][0];
    t4nyz = c1*LC->axis_lc[8][0];

    c.gmt1x = t1x + s1 * I_t * (kvec[0]*LC->axis_gc_x[0][3] + kvec[1]*LC->axis_gc_y[0][3] + kvec[2]*LC->axis_gc_z[0][3]);
    c.gmt2x = t2x + s1 * I_t * (kvec[0]*LC->axis_gc_x[0][2] + kvec[1]*LC->axis_gc_y[0][2] + kvec[2]*LC->axis_gc_z[0][2]);
    c.gmt3x = t3x + s1 * I_t * (kvec[0]*LC->axis_gc_x[0][1] + kvec[1]*LC->axis_gc_y[0][1] + kvec[2]*LC->axis_gc_z[0][1]);
    c.gmt4x = t4x + s1 * I_t * (kvec[0]*LC->axis_gc_x[0][0] + kvec[1]*LC->axis_gc_y[0][0] + kvec[2]*LC->axis_gc_z[0][0]);

    c.gpt1x = t1x - s1 * I_t * (kvec[0]*LC->axis_gc_x[0][3] + kvec[1]*LC->axis_gc_y[0][3] + kvec[2]*LC->axis_gc_z[0][3]);
    c.gpt2x = t2x - s1 * I_t * (kvec[0]*LC->axis_gc_x[0][2] + kvec[1]*LC->axis_gc_y[0][2] + kvec[2]*LC->axis_gc_z[0][2]);
    c.gpt3x = t3x - s1 * I_t * (kvec[0]*LC->axis_gc_x[0][1] + kvec[1]*LC->axis_gc_y[0][1] + kvec[2]*LC->axis_gc_z[0][1]);
    c.gpt4x = t4x - s1 * I_t * (kvec[0]*LC->axis_gc_x[0][0] + kvec[1]*LC->axis_gc_y[0][0] + kvec[2]*LC->axis_gc_z[0][0]);

    c.gmt1y = t1y + s1 * I_t * (kvec[0]*LC->axis_gc_x[1][3] + kvec[1]*LC->axis_gc_y[1][3] + kvec[2]*LC->axis_gc_z[1][3]);
    c.gmt2y = t2y + s1 * I_t * (kvec[0]*LC->axis_gc_x[1][2] + kvec[1]*LC->axis_gc_y[1][2] + kvec[2]*LC->axis_gc_z[1][2]);
    c.gmt3y = t3y + s1 * I_t * (kvec[0]*LC->axis_gc_x[1][1] + kvec[1]*LC->axis_gc_y[1][1] + kvec[2]*LC->axis_gc_z[1][1]);
    c.gmt4y = t4y + s1 * I_t * (kvec[0]*LC->axis_gc_x[1][0] + kvec[1]*LC->axis_gc_y[1][0] + kvec[2]*LC->axis_gc_z[1][0]);

    c.gpt1y = t1y - s1 * I_t * (kvec[0]*LC->axis_gc_x[1][3] + kvec[1]*LC->axis_gc_y[1][3] + kvec[2]*LC->axis_gc_z[1][3]);
    c.gpt2y = t2y - s1 * I_t * (kvec[0]*LC->axis_gc_x[1][2] + kvec[1]*LC->axis_gc_y[1][2] + kvec[2]*LC->axis_gc_z[1][2]);
    c.gpt3y = t3y - s1 * I_t * (kvec[0]*LC->axis_gc_x[1][1] + kvec[1]*LC->axis_gc_y[1][1] + kvec[2]*LC->axis_gc_z[1][1]);
    c.gpt4y = t4y - s1 * I_t * (kvec[0]*LC->axis_gc_x[1][0] + kvec[1]*LC->axis_gc_y[1][0] + kvec[2]*LC->axis_gc_z[1][0]);

    c.gmt1z = t1z + s1 * I_t * (kvec[0]*LC->axis_gc_x[2][3] + kvec[1]*LC->axis_gc_y[2][3] + kvec[2]*LC->axis_gc_z[2][3]);
    c.gmt2z = t2z + s1 * I_t * (kvec[0]*LC->axis_gc_x[2][2] + kvec[1]*LC->axis_gc_y[2][2] + kvec[2]*LC->axis_gc_z[2][2]);
    c.gmt3z = t3z + s1 * I_t * (kvec[0]*LC->axis_gc_x[2][1] + kvec[1]*LC->axis_gc_y[2][1] + kvec[2]*LC->axis_gc_z[2][1]);
    c.gmt4z = t4z + s1 * I_t * (kvec[0]*LC->axis_gc_x[2][0] + kvec[1]*LC->axis_gc_y[2][0] + kvec[2]*LC->axis_gc_z[2][0]);

    c.gpt1z = t1z - s1 * I_t * (kvec[0]*LC->axis_gc_x[2][3] + kvec[1]*LC->axis_gc_y[2][3] + kvec[2]*LC->axis_gc_z[2][3]);
    c.gpt2z = t2z - s1 * I_t * (kvec[0]*LC->axis_gc_x[2][2] + kvec[1]*LC->axis_gc_y[2][2] + kvec[2]*LC->axis_gc_z[2][2]);
    c.gpt3z = t3z - s1 * I_t * (kvec[0]*LC->axis_gc_x[2][1] + kvec[1]*LC->axis_gc_y[2][1] + kvec[2]*LC->axis_gc_z[2][1]);
    c.gpt4z = t4z - s1 * I_t * (kvec[0]*LC->axis_gc_x[2][0] + kvec[1]*LC->axis_gc_y[2][0] + kvec[2]*LC->axis_gc_z[2][0]);

    c.gmt1xy = t1xy + s1 * I_t * (kvec[0]*LC->axis_gc_x[3][3] + kvec[1]*LC->axis_gc_y[3][3] + kvec[2]*LC->axis_gc_z[3][3]);
    c.gmt2xy = t2xy + s1 * I_t * (kvec[0]*LC->axis_gc_x[3][2] + kvec[1]*LC->axis_gc_y[3][2] + kvec[2]*LC->axis_gc_z[3][2]);
    c.gmt3xy = t3xy + s1 * I_t * (kvec[0]*LC->axis_gc_x[3][1] + kvec[1]*LC->axis_gc_y[3][1] + kvec[2]*LC->axis_gc_z[3][1]);
    c.gmt4xy = t4xy + s1 * I_t * (kvec[0]*LC->axis_gc_x[3][0] + kvec[1]*LC->axis_gc_y[3][0] + kvec[2]*LC->axis_gc_z[3][0]);

    c.gpt1xy = t1xy - s1 * I_t * (kvec[0]*LC->axis_gc_x[3][3] + kvec[1]*LC->axis_gc_y[3][3] + kvec[2]*LC->axis_gc_z[3][3]);
    c.gpt2xy = t2xy - s1 * I_t * (kvec[0]*LC->axis_gc_x[3][2] + kvec[1]*LC->axis_gc_y[3][2] + kvec[2]*LC->axis_gc_z[3][2]);
    c.gpt3xy = t3xy - s1 * I_t * (kvec[0]*LC->axis_gc_x[3][1] + kvec[1]*LC->axis_gc_y[3][1] + kvec[2]*LC->axis_gc_z[3][1]);
    c.gpt4xy = t4xy - s1 * I_t * (kvec[0]*LC->axis_gc_x[3][0] + kvec[1]*LC->axis_gc_y[3][0] + kvec[2]*LC->axis_gc_z[3][0]);

    c.gmt1nxy = t1nxy + s1 * I_t * (kvec[0]*LC->axis_gc_x[6][3] + kvec[1]*LC->axis_gc_y[6][3] + kvec[2]*LC->axis_gc_z[6][3]);
    c.gmt2nxy = t2nxy + s1 * I_t * (kvec[0]*LC->axis_gc_x[6][2] + kvec[1]*LC->axis_gc_y[6][2] + kvec[2]*LC->axis_gc_z[6][2]);
    c.gmt3nxy = t3nxy + s1 * I_t * (kvec[0]*LC->axis_gc_x[6][1] + kvec[1]*LC->axis_gc_y[6][1] + kvec[2]*LC->axis_gc_z[6][1]);
    c.gmt4nxy = t4nxy + s1 * I_t * (kvec[0]*LC->axis_gc_x[6][0] + kvec[1]*LC->axis_gc_y[6][0] + kvec[2]*LC->axis_gc_z[6][0]);

    c.gpt1nxy = t1nxy - s1 * I_t * (kvec[0]*LC->axis_gc_x[6][3] + kvec[1]*LC->axis_gc_y[6][3] + kvec[2]*LC->axis_gc_z[6][3]);
    c.gpt2nxy = t2nxy - s1 * I_t * (kvec[0]*LC->axis_gc_x[6][2] + kvec[1]*LC->axis_gc_y[6][2] + kvec[2]*LC->axis_gc_z[6][2]);
    c.gpt3nxy = t3nxy - s1 * I_t * (kvec[0]*LC->axis_gc_x[6][1] + kvec[1]*LC->axis_gc_y[6][1] + kvec[2]*LC->axis_gc_z[6][1]);
    c.gpt4nxy = t4nxy - s1 * I_t * (kvec[0]*LC->axis_gc_x[6][0] + kvec[1]*LC->axis_gc_y[6][0] + kvec[2]*LC->axis_gc_z[6][0]);

    c.gmt1xz = t1xz + s1 * I_t * (kvec[0]*LC->axis_gc_x[4][3] + kvec[1]*LC->axis_gc_y[4][3] + kvec[2]*LC->axis_gc_z[4][3]);
    c.gmt2xz = t2xz + s1 * I_t * (kvec[0]*LC->axis_gc_x[4][2] + kvec[1]*LC->axis_gc_y[4][2] + kvec[2]*LC->axis_gc_z[4][2]);
    c.gmt3xz = t3xz + s1 * I_t * (kvec[0]*LC->axis_gc_x[4][1] + kvec[1]*LC->axis_gc_y[4][1] + kvec[2]*LC->axis_gc_z[4][1]);
    c.gmt4xz = t4xz + s1 * I_t * (kvec[0]*LC->axis_gc_x[4][0] + kvec[1]*LC->axis_gc_y[4][0] + kvec[2]*LC->axis_gc_z[4][0]);

    c.gpt1xz = t1xz - s1 * I_t * (kvec[0]*LC->axis_gc_x[4][3] + kvec[1]*LC->axis_gc_y[4][3] + kvec[2]*LC->axis_gc_z[4][3]);
    c.gpt2xz = t2xz - s1 * I_t * (kvec[0]*LC->axis_gc_x[4][2] + kvec[1]*LC->axis_gc_y[4][2] + kvec[2]*LC->axis_gc_z[4][2]);
    c.gpt3xz = t3xz - s1 * I_t * (kvec[0]*LC->axis_gc_x[4][1] + kvec[1]*LC->axis_gc_y[4][1] + kvec[2]*LC->axis_gc_z[4][1]);
    c.gpt4xz = t4xz - s1 * I_t * (kvec[0]*LC->axis_gc_x[4][0] + kvec[1]*LC->axis_gc_y[4][0] + kvec[2]*LC->axis_gc_z[4][0]);

    c.gmt1nxz = t1nxz + s1 * I_t * (kvec[0]*LC->axis_gc_x[7][3] + kvec[1]*LC->axis_gc_y[7][3] + kvec[2]*LC->axis_gc_z[7][3]);
    c.gmt2nxz = t2nxz + s1 * I_t * (kvec[0]*LC->axis_gc_x[7][2] + kvec[1]*LC->axis_gc_y[7][2] + kvec[2]*LC->axis_gc_z[7][2]);
    c.gmt3nxz = t3nxz + s1 * I_t * (kvec[0]*LC->axis_gc_x[7][1] + kvec[1]*LC->axis_gc_y[7][1] + kvec[2]*LC->axis_gc_z[7][1]);
    c.gmt4nxz = t4nxz + s1 * I_t * (kvec[0]*LC->axis_gc_x[7][0] + kvec[1]*LC->axis_gc_y[7][0] + kvec[2]*LC->axis_gc_z[7][0]);

    c.gpt1nxz = t1nxz - s1 * I_t * (kvec[0]*LC->axis_gc_x[7][3] + kvec[1]*LC->axis_gc_y[7][3] + kvec[2]*LC->axis_gc_z[7][3]);
    c.gpt2nxz = t2nxz - s1 * I_t * (kvec[0]*LC->axis_gc_x[7][2] + kvec[1]*LC->axis_gc_y[7][2] + kvec[2]*LC->axis_gc_z[7][2]);
    c.gpt3nxz = t3nxz - s1 * I_t * (kvec[0]*LC->axis_gc_x[7][1] + kvec[1]*LC->axis_gc_y[7][1] + kvec[2]*LC->axis_gc_z[7][1]);
    c.gpt4nxz = t4nxz - s1 * I_t * (kvec[0]*LC->axis_gc_x[7][0] + kvec[1]*LC->axis_gc_y[7][0] + kvec[2]*LC->axis_gc_z[7][0]);

    c.gmt1yz = t1yz + s1 * I_t * (kvec[0]*LC->axis_gc_x[5][3] + kvec[1]*LC->axis_gc_y[5][3] + kvec[2]*LC->axis_gc_z[5][3]);
    c.gmt2yz = t2yz + s1 * I_t * (kvec[0]*LC->axis_gc_x[5][2] + kvec[1]*LC->axis_gc_y[5][2] + kvec[2]*LC->axis_gc_z[5][2]);
    c.gmt3yz = t3yz + s1 * I_t * (kvec[0]*LC->axis_gc_x[5][1] + kvec[1]*LC->axis_gc_y[5][1] + kvec[2]*LC->axis_gc_z[5][1]);
    c.gmt4yz = t4yz + s1 * I_t * (kvec[0]*LC->axis_gc_x[5][0] + kvec[1]*LC->axis_gc_y[5][0] + kvec[2]*LC->axis_gc_z[5][0]);

    c.gpt1yz = t1yz - s1 * I_t * (kvec[0]*LC->axis_gc_x[5][3] + kvec[1]*LC->axis_gc_y[5][3] + kvec[2]*LC->axis_gc_z[5][3]);
    c.gpt2yz = t2yz - s1 * I_t * (kvec[0]*LC->axis_gc_x[5][2] + kvec[1]*LC->axis_gc_y[5][2] + kvec[2]*LC->axis_gc_z[5][2]);
    c.gpt3yz = t3yz - s1 * I_t * (kvec[0]*LC->axis_gc_x[5][1] + kvec[1]*LC->axis_gc_y[5][1] + kvec[2]*LC->axis_gc_z[5][1]);
    c.gpt4yz = t4yz - s1 * I_t * (kvec[0]*LC->axis_gc_x[5][0] + kvec[1]*LC->axis_gc_y[5][0] + kvec[2]*LC->axis_gc_z[5][0]);

    c.gmt1nyz = t1nyz + s1 * I_t * (kvec[0]*LC->axis_gc_x[8][3] + kvec[1]*LC->axis_gc_y[8][3] + kvec[2]*LC->axis_gc_z[8][3]);
    c.gmt2nyz = t2nyz + s1 * I_t * (kvec[0]*LC->axis_gc_x[8][2] + kvec[1]*LC->axis_gc_y[8][2] + kvec[2]*LC->axis_gc_z[8][2]);
    c.gmt3nyz = t3nyz + s1 * I_t * (kvec[0]*LC->axis_gc_x[8][1] + kvec[1]*LC->axis_gc_y[8][1] + kvec[2]*LC->axis_gc_z[8][1]);
    c.gmt4nyz = t4nyz + s1 * I_t * (kvec[0]*LC->axis_gc_x[8][0] + kvec[1]*LC->axis_gc_y[8][0] + kvec[2]*LC->axis_gc_z[8][0]);

    c.gpt1nyz = t1nyz - s1 * I_t * (kvec[0]*LC->axis_gc_x[8][3] + kvec[1]*LC->axis_gc_y[8][3] + kvec[2]*LC->axis_gc_z[8][3]);
    c.gpt2nyz = t2nyz - s1 * I_t * (kvec[0]*LC->axis_gc_x[8][2] + kvec[1]*LC->axis_gc_y[8][2] + kvec[2]*LC->axis_gc_z[8][2]);
    c.gpt3nyz = t3nyz - s1 * I_t * (kvec[0]*LC->axis_gc_x[8][1] + kvec[1]*LC->axis_gc_y[8][1] + kvec[2]*LC->axis_gc_z[8][1]);
    c.gpt4nyz = t4nyz - s1 * I_t * (kvec[0]*LC->axis_gc_x[8][0] + kvec[1]*LC->axis_gc_y[8][0] + kvec[2]*LC->axis_gc_z[8][0]);

    // NULL b means we just want the diagonal component.
    if(b == NULL) return (double)std::real(th2);

    // x,y, and z axes are required by all lattice types
    for (int ix = 4; ix < dimx + 4; ix++)
    {

        for (int iy = 4; iy < dimy + 4; iy++)
        {

            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
            // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
            // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] = th2 * A[iz] +
                    c.gpt1z * A[iz + 1] + c.gmt1z * A[iz - 1] +
                    c.gpt2z * A[iz + 2] + c.gmt2z * A[iz - 2] +
                    c.gpt3z * A[iz + 3] + c.gmt3z * A[iz - 3] +
                    c.gpt4z * A[iz + 4] + c.gmt4z * A[iz - 4];
            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    c.gpt1y * A[iz + iys] + c.gmt1y * A[iz - iys] +
                    c.gpt2y * A[iz + 2*iys] + c.gmt2y * A[iz - 2*iys] +
                    c.gpt3y * A[iz + 3*iys] + c.gmt3y * A[iz - 3*iys] +
                    c.gpt4y * A[iz + 4*iys] + c.gmt4y * A[iz - 4*iys];

            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    c.gpt1x * A[iz + ixs] + c.gmt1x * A[iz - ixs] +
                    c.gpt2x * A[iz + 2*ixs] + c.gmt2x * A[iz - 2*ixs] +
                    c.gpt3x * A[iz + 3*ixs] + c.gmt3x * A[iz - 3*ixs] +
                    c.gpt4x * A[iz + 4*ixs] + c.gmt4x * A[iz - 4*ixs];
            }                   /* end for */
        }
    }

    /* Quick return for orthogonal axis cases */
    if(ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == CUBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
        return (double)std::real(th2);

    for (int ix = 4; ix < dimx + 4; ix++)
    {

        for (int iy = 4; iy < dimy + 4; iy++)
        {

            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
            if(LC->include_axis[3])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1xy * A[iz + ixs + iys] + c.gmt1xy * A[iz - ixs - iys] +
                        c.gpt2xy * A[iz + 2*ixs + 2*iys] + c.gmt2xy * A[iz - 2*ixs - 2*iys] +
                        c.gpt3xy * A[iz + 3*ixs + 3*iys] + c.gmt3xy * A[iz - 3*ixs - 3*iys] +
                        c.gpt4xy * A[iz + 4*ixs + 4*iys] + c.gmt4xy * A[iz - 4*ixs - 4*iys];
                }                   /* end for */
            }
            if(LC->include_axis[4])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1xz * A[iz + ixs + 1] + c.gmt1xz * A[iz - ixs - 1] +
                        c.gpt2xz * A[iz + 2*ixs + 2] + c.gmt2xz * A[iz - 2*ixs - 2] +
                        c.gpt3xz * A[iz + 3*ixs + 3] + c.gmt3xz * A[iz - 3*ixs - 3] +
                        c.gpt4xz * A[iz + 4*ixs + 4] + c.gmt4xz * A[iz - 4*ixs - 4];
                }                   /* end for */
            }
            if(LC->include_axis[5])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1yz * A[iz + iys + 1] + c.gmt1yz * A[iz - iys - 1] +
                        c.gpt2yz * A[iz + 2*iys + 2] + c.gmt2yz * A[iz - 2*iys - 2] +
                        c.gpt3yz * A[iz + 3*iys + 3] + c.gmt3yz * A[iz - 3*iys - 3] +
                        c.gpt4yz * A[iz + 4*iys + 4] + c.gmt4yz * A[iz - 4*iys - 4];
                }                   /* end for */
            }
            if(LC->include_axis[6])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1nxy * A[iz - ixs + iys] + c.gmt1nxy * A[iz + ixs - iys] +
                        c.gpt2nxy * A[iz - 2*ixs + 2*iys] + c.gmt2nxy * A[iz + 2*ixs - 2*iys] +
                        c.gpt3nxy * A[iz - 3*ixs + 3*iys] + c.gmt3nxy * A[iz + 3*ixs - 3*iys] +
                        c.gpt4nxy * A[iz - 4*ixs + 4*iys] + c.gmt4nxy * A[iz + 4*ixs - 4*iys];
                }                   /* end for */
            }
            if(LC->include_axis[7])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1nxz * A[iz - ixs + 1] + c.gmt1nxz * A[iz + ixs - 1] +
                        c.gpt2nxz * A[iz - 2*ixs + 2] + c.gmt2nxz * A[iz + 2*ixs - 2] +
                        c.gpt3nxz * A[iz - 3*ixs + 3] + c.gmt3nxz * A[iz + 3*ixs - 3] +
                        c.gpt4nxz * A[iz - 4*ixs + 4] + c.gmt4nxz * A[iz + 4*ixs - 4];
                }                   /* end for */
            }
            if(LC->include_axis[8])
            {
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        c.gpt1nyz * A[iz - iys + 1] + c.gmt1nyz * A[iz + iys - 1] +
                        c.gpt2nyz * A[iz - 2*iys + 2] + c.gmt2nyz * A[iz + 2*iys - 2] +
                        c.gpt3nyz * A[iz - 3*iys + 3] + c.gmt3nyz * A[iz + 3*iys - 3] +
                        c.gpt4nyz * A[iz - 4*iys + 4] + c.gmt4nyz * A[iz + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    /* Return the diagonal component of the operator */
    return (double)std::real(th2);


} /* end app8_combined_general */
