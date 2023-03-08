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


template double FiniteDiff::app10_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app10_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);

template void FiniteDiff::app_gradient_tenth<float> (float *, float *, float *, float *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_tenth<double> (double *, double *, double *, double *, int, int, int, double, double, double);
template void FiniteDiff::app_gradient_tenth<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int, double , double , double );
template void FiniteDiff::app_gradient_tenth<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int, double , double , double );


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



#include "rmg_complex.h"


template double FiniteDiff::app12_del2<float>(float *, float *, int, int, int, double, double, double);
template double FiniteDiff::app12_del2<double>(double *, double *, int, int, int, double, double, double);
template double FiniteDiff::app12_del2<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double);
template double FiniteDiff::app12_del2<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double);


template <typename RmgType>
double FiniteDiff::app12_del2(RmgType * a, RmgType * b, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{

    int ibrav = L->get_ibrav_type();
    int iz, ix, iy;
    int ixs, iys, ix1, iy1;

    ixs = (dimy + 12) * (dimz + 12);
    iys = (dimz + 12);
    ix1 = dimy * dimz;
    iy1 = dimz;

    // nine and seven point stencils, 2nd derivative, extrapolated
    int ic = 6;
    double x[14], w1[14], w2[14];
    for(int i=0;i<14;i++) x[i] = (double)i;
    gen_weights(13, 2, (double)ic, x, w1);
    gen_weights(11, 2, (double)(ic-1), x, w2);
    double hf = 1.0, c1, c2=0.0;
    if(ibrav == HEXAGONAL || ibrav == HEXAGONAL2) hf = 2.0/3.0;

    double h2x = gridhx * gridhx * L->get_xside() * L->get_xside();
    double h2y = gridhy * gridhy * L->get_yside() * L->get_yside();
    double h2z = gridhz * gridhz * L->get_zside() * L->get_zside();

    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    double th2 = (c1*w1[ic] - c2*w2[ic-1]) / h2x;
    RmgType t1x ((c1*w1[ic+1] - c2*w2[ic]) * hf / h2x);
    RmgType t2x ((c1*w1[ic+2] - c2*w2[ic+1]) * hf / h2x);
    RmgType t3x ((c1*w1[ic+3] - c2*w2[ic+2]) * hf / h2x);
    RmgType t4x ((c1*w1[ic+4] - c2*w2[ic+3]) * hf / h2x);
    RmgType t5x ((c1*w1[ic+5] - c2*w2[ic+4]) * hf / h2x);
    RmgType t6x (c1*w1[ic+6] * hf / h2x);

    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) /  h2y;
    RmgType t1y ((c1*w1[ic+1] - c2*w2[ic]) / h2y);
    RmgType t2y ((c1*w1[ic+2] - c2*w2[ic+1]) / h2y);
    RmgType t3y ((c1*w1[ic+3] - c2*w2[ic+2]) / h2y);
    RmgType t4y ((c1*w1[ic+4] - c2*w2[ic+3]) / h2y);
    RmgType t5y ((c1*w1[ic+5] - c2*w2[ic+4]) / h2y);
    RmgType t6y (c1*w1[ic+6] / h2y);

    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    th2 += (c1*w1[ic] - c2*w2[ic-1]) /  h2z;
    RmgType t1z ((c1*w1[ic+1] - c2*w2[ic]) / h2z);
    RmgType t2z ((c1*w1[ic+2] - c2*w2[ic+1]) / h2z);
    RmgType t3z ((c1*w1[ic+3] - c2*w2[ic+2]) / h2z);
    RmgType t4z ((c1*w1[ic+4] - c2*w2[ic+3]) / h2z);
    RmgType t5z ((c1*w1[ic+5] - c2*w2[ic+4]) / h2z);
    RmgType t6z (c1*w1[ic+6] / h2z);
    RmgType t0 (th2);

    // NULL b means we just want the diagonal component.
    if(b == NULL) return (double)std::real(t0);

    switch(ibrav)
    {
        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:
        case TETRAGONAL_PRIMITIVE:

            for (ix = 6; ix < dimx + 6; ix++)
            {

                for (iy = 6; iy < dimy + 6; iy++)
                {

                    for (iz = 6; iz < dimz + 6; iz++)
                    {

                        b[(ix - 6) * ix1 + (iy - 6) * iy1 + iz - 6] =
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
                            t6x * a[(ix - 6) * ixs + iy * iys + iz] +
                            t6x * a[(ix + 6) * ixs + iy * iys + iz] +

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
                            t6y * a[ix * ixs + (iy - 6) * iys + iz] +
                            t6y * a[ix * ixs + (iy + 6) * iys + iz] +

                            t1z * a[ix * ixs + iy * iys + iz - 1] +
                            t1z * a[ix * ixs + iy * iys + iz + 1] +
                            t2z * a[ix * ixs + iy * iys + iz - 2] +
                            t2z * a[ix * ixs + iy * iys + iz + 2] +
                            t3z * a[ix * ixs + iy * iys + iz - 3] +
                            t3z * a[ix * ixs + iy * iys + iz + 3] +
                            t4z * a[ix * ixs + iy * iys + iz - 4] +
                            t4z * a[ix * ixs + iy * iys + iz + 4] +
                            t5z * a[ix * ixs + iy * iys + iz - 5] +
                            t5z * a[ix * ixs + iy * iys + iz + 5] +
                            t6z * a[ix * ixs + iy * iys + iz - 6] +
                            t6z * a[ix * ixs + iy * iys + iz + 6];

                    }                   /* end for */
                }                       /* end for */
            }                           /* end for */
            break;

        case HEXAGONAL2:
        case HEXAGONAL:

            for (int ix = 6; ix < dimx + 6; ix++)
            {

                for (int iy = 6; iy < dimy + 6; iy++)
                {

                    RmgType *A = &a[iy*iys + ix*ixs];
                    RmgType *B = &b[(iy - 6)*dimz + (ix - 6)*dimy*dimz - 6];
                    // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
                    for (int iz = 6; iz < dimz + 6; iz++)
                    {
                        B[iz] = t0 * A[iz] +
                                t1z * (A[iz + 1] + A[iz - 1]) +
                                t2z * (A[iz + 2] + A[iz - 2]) +
                                t3z * (A[iz + 3] + A[iz - 3]) +
                                t4z * (A[iz + 4] + A[iz - 4]) +
                                t5z * (A[iz + 5] + A[iz - 5]) +
                                t6z * (A[iz + 6] + A[iz - 6]);
                    }

                    for (int iz = 6; iz < dimz + 6; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + iys] + A[iz - iys]) +
                                t2x * (A[iz + 2*iys] + A[iz - 2*iys]) +
                                t3x * (A[iz + 3*iys] + A[iz - 3*iys]) +
                                t4x * (A[iz + 4*iys] + A[iz - 4*iys]) +
                                t5x * (A[iz + 5*iys] + A[iz - 5*iys]) +
                                t6x * (A[iz + 6*iys] + A[iz - 6*iys]);
                    }

                    for (int iz = 6; iz < dimz + 6; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + ixs] + A[iz - ixs]) +
                                t2x * (A[iz + 2*ixs] + A[iz - 2*ixs]) +
                                t3x * (A[iz + 3*ixs] + A[iz - 3*ixs]) +
                                t4x * (A[iz + 4*ixs] + A[iz - 4*ixs]) +
                                t5x * (A[iz + 5*ixs] + A[iz - 5*ixs]) +
                                t6x * (A[iz + 6*ixs] + A[iz - 6*ixs]);
                    }                   /* end for */

                    for (int iz = 6; iz < dimz + 6; iz++)
                    {
                        B[iz] +=
                                t1x * (A[iz + ixs + iys] + A[iz - ixs - iys]) +
                                t2x * (A[iz + 2*ixs + 2*iys] + A[iz - 2*ixs - 2*iys]) +
                                t3x * (A[iz + 3*ixs + 3*iys] + A[iz - 3*ixs - 3*iys]) +
                                t4x * (A[iz + 4*ixs + 4*iys] + A[iz - 4*ixs - 4*iys]) +
                                t5x * (A[iz + 5*ixs + 5*iys] + A[iz - 5*ixs - 5*iys]) +
                                t6x * (A[iz + 6*ixs + 6*iys] + A[iz - 6*ixs - 6*iys]);
                    }                   /* end for */

                }
            }

            break;

        default:
            rmg_error_handler (__FILE__, __LINE__, "Lattice type not implemented");

    }


    /* Return the diagonal component of the operator */
    return (double)std::real(t0);

} /* end app12_del2 */

