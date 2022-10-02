/*
 *
 * Copyright (c) 1995,2011,2014,2022 Emil Briggs
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
#include "rmg_complex.h"


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


template void FiniteDiff::fd_gradient_coeffs<float>(int , double, int , float *, float *, float *);
template void FiniteDiff::fd_gradient_coeffs<double>(int , double, int , double *, double *, double *);
template void FiniteDiff::fd_gradient_coeffs<std::complex<double>>(int , double, int , std::complex<double> *, std::complex<double> *, std::complex<double> *);

template void FiniteDiff::fd_combined_coeffs<float>(int, double, int, float *, float *, double *);
template void FiniteDiff::fd_combined_coeffs<double>(int, double, int, double *, double *, double *);
template void FiniteDiff::fd_combined_coeffs<std::complex<float>>(int, double, int, std::complex<float> *, std::complex<float> *, double *);
template void FiniteDiff::fd_combined_coeffs<std::complex<double>>(int, double, int, std::complex<double> *, std::complex<double> *, double *);


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
    RmgType cm[8], cp[8];
    RmgType cpx[8], cmx[8], cpy[8], cmy[8], cpz[8], cmz[8];
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
    RmgType cxx[8], cxy[8], cxz[8];
    RmgType cyx[8], cyy[8], cyz[8];
    RmgType czx[8], czy[8], czz[8];
    RmgType cx[8], cy[8], cz[8];
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

                bgy[iz] = -czy[0] * A[iz + 1] + czy[0] * A[iz - 1];
                if constexpr(order >= 4)
                    bgy[iz] += -czy[1] * A[iz + 2] + czy[1] * A[iz - 2];
                if constexpr(order >= 6)
                    bgy[iz] += -czy[2] * A[iz + 3] + czy[2] * A[iz - 3];
                if constexpr(order >= 8)
                    bgy[iz] += -czy[3] * A[iz + 4] + czy[3] * A[iz - 4];

                bgz[iz] = -czz[0] * A[iz + 1] + czz[0] * A[iz - 1];
                if constexpr(order >= 4)
                    bgz[iz] += -czz[1] * A[iz + 2] + czz[1] * A[iz - 2];
                if constexpr(order >= 6)
                    bgz[iz] += -czz[2] * A[iz + 3] + czz[2] * A[iz - 3];
                if constexpr(order >= 8)
                    bgz[iz] += -czz[3] * A[iz + 4] + czz[3] * A[iz - 4];

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

                bgy[iz] += -cyy[0] * A[iz + iys] + cyy[0] * A[iz - iys];
                if constexpr(order >= 4)
                    bgy[iz] += -cyy[1] * A[iz + 2*iys] + cyy[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    bgy[iz] += -cyy[2] * A[iz + 3*iys] + cyy[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    bgy[iz] += -cyy[3] * A[iz + 4*iys] + cyy[3] * A[iz - 4*iys];

                bgz[iz] += -cyz[0] * A[iz + iys] + cyz[0] * A[iz - iys];
                if constexpr(order >= 4)
                    bgz[iz] += -cyz[1] * A[iz + 2*iys] + cyz[1] * A[iz - 2*iys];
                if constexpr(order >= 6)
                    bgz[iz] += -cyz[2] * A[iz + 3*iys] + cyz[2] * A[iz - 3*iys];
                if constexpr(order >= 8)
                    bgz[iz] += -cyz[3] * A[iz + 4*iys] + cyz[3] * A[iz - 4*iys];
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

                bgy[iz] += -cxy[0] * A[iz + ixs] + cxy[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    bgy[iz] +=-cxy[1] * A[iz + 2*ixs] + cxy[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    bgy[iz] +=-cxy[2] * A[iz + 3*ixs] + cxy[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    bgy[iz] +=-cxy[3] * A[iz + 4*ixs] + cxy[3] * A[iz - 4*ixs];

                bgz[iz] += -cxz[0] * A[iz + ixs] + cxz[0] * A[iz - ixs];
                if constexpr(order >= 4)
                    bgz[iz] += -cxz[1] * A[iz + 2*ixs] + cxz[1] * A[iz - 2*ixs];
                if constexpr(order >= 6)
                    bgz[iz] += -cxz[2] * A[iz + 3*ixs] + cxz[2] * A[iz - 3*ixs];
                if constexpr(order >= 8)
                    bgz[iz] += -cxz[3] * A[iz + 4*ixs] + cxz[3] * A[iz - 4*ixs];

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

                    bgy[iz] += -cy[0] * A[iz + ixs + iys] + cy[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2*iys] + cy[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3*iys] + cy[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4*iys] + cy[3] * A[iz - 4*ixs - 4*iys];

                    bgx[iz] += -cx[0] * A[iz + ixs + iys] + cx[0] * A[iz - ixs - iys];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2*iys] + cx[1] * A[iz - 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3*iys] + cx[2] * A[iz - 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4*iys] + cx[3] * A[iz - 4*ixs - 4*iys];
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

                    bgy[iz] += -cy[0] * A[iz + ixs + 1] + cy[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2] + cy[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3] + cy[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4] + cy[3] * A[iz - 4*ixs - 4];

                    bgx[iz] += -cx[0] * A[iz + ixs + 1] + cx[0] * A[iz - ixs - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2] + cx[1] * A[iz - 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3] + cx[2] * A[iz - 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4] + cx[3] * A[iz - 4*ixs - 4];
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

                    bgy[iz] += -cy[0] * A[iz + iys + 1] + cy[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*iys + 2] + cy[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*iys + 3] + cy[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*iys + 4] + cy[3] * A[iz - 4*iys - 4];

                    bgx[iz] += -cx[0] * A[iz + iys + 1] + cx[0] * A[iz - iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*iys + 2] + cx[1] * A[iz - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*iys + 3] + cx[2] * A[iz - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*iys + 4] + cx[3] * A[iz - 4*iys - 4];
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

                    bgy[iz] += -cy[0] * A[iz - ixs + iys] + cy[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs + 2*iys] + cy[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs + 3*iys] + cy[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs + 4*iys] + cy[3] * A[iz + 4*ixs - 4*iys];

                    bgx[iz] += -cx[0] * A[iz - ixs + iys] + cx[0] * A[iz + ixs - iys];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs + 2*iys] + cx[1] * A[iz + 2*ixs - 2*iys];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs + 3*iys] + cx[2] * A[iz + 3*ixs - 3*iys];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs + 4*iys] + cx[3] * A[iz + 4*ixs - 4*iys];
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

                    bgy[iz] += -cy[0] * A[iz - ixs + 1] + cy[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs + 2] + cy[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs + 3] + cy[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs + 4] + cy[3] * A[iz + 4*ixs - 4];

                    bgx[iz] += -cx[0] * A[iz - ixs + 1] + cx[0] * A[iz + ixs - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs + 2] + cx[1] * A[iz + 2*ixs - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs + 3] + cx[2] * A[iz + 3*ixs - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs + 4] + cx[3] * A[iz + 4*ixs - 4];
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

                    bgy[iz] += -cy[0] * A[iz - iys + 1] + cy[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*iys + 2] + cy[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*iys + 3] + cy[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*iys + 4] + cy[3] * A[iz + 4*iys - 4];

                    bgx[iz] += -cx[0] * A[iz - iys + 1] + cx[0] * A[iz + iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*iys + 2] + cx[1] * A[iz + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*iys + 3] + cx[2] * A[iz + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*iys + 4] + cx[3] * A[iz + 4*iys - 4];
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

                    bgy[iz] += -cy[0] * A[iz + 1*ixs + 1*iys + 1] + cy[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs + 2*iys + 2] + cy[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs + 3*iys + 3] + cy[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs + 4*iys + 4] + cy[3] * A[iz - 4*ixs - 4*iys - 4];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs + 1*iys + 1] + cx[0] * A[iz - 1*ixs - 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs + 2*iys + 2] + cx[1] * A[iz - 2*ixs - 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs + 3*iys + 3] + cx[2] * A[iz - 3*ixs - 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs + 4*iys + 4] + cx[3] * A[iz - 4*ixs - 4*iys - 4];
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

                    bgy[iz] += -cy[0] * A[iz - 1*ixs - 1*iys + 1] + cy[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz - 2*ixs - 2*iys + 2] + cy[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz - 3*ixs - 3*iys + 3] + cy[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz - 4*ixs - 4*iys + 4] + cy[3] * A[iz + 4*ixs + 4*iys - 4];

                    bgx[iz] += -cx[0] * A[iz - 1*ixs - 1*iys + 1] + cx[0] * A[iz + 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz - 2*ixs - 2*iys + 2] + cx[1] * A[iz + 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz - 3*ixs - 3*iys + 3] + cx[2] * A[iz + 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz - 4*ixs - 4*iys + 4] + cx[3] * A[iz + 4*ixs + 4*iys - 4];
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

                    bgy[iz] += -cy[0] * A[iz + 1*ixs - 1*iys + 1] + cy[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs - 2*iys + 2] + cy[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs - 3*iys + 3] + cy[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs - 4*iys + 4] + cy[3] * A[iz - 4*ixs + 4*iys - 4];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs - 1*iys + 1] + cx[0] * A[iz - 1*ixs + 1*iys - 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs - 2*iys + 2] + cx[1] * A[iz - 2*ixs + 2*iys - 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs - 3*iys + 3] + cx[2] * A[iz - 3*ixs + 3*iys - 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs - 4*iys + 4] + cx[3] * A[iz - 4*ixs + 4*iys - 4];
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

                    bgy[iz] += -cy[0] * A[iz + 1*ixs - 1*iys - 1] + cy[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        bgy[iz] += -cy[1] * A[iz + 2*ixs - 2*iys - 2] + cy[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        bgy[iz] += -cy[2] * A[iz + 3*ixs - 3*iys - 3] + cy[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        bgy[iz] += -cy[3] * A[iz + 4*ixs - 4*iys - 4] + cy[3] * A[iz - 4*ixs + 4*iys + 4];

                    bgx[iz] += -cx[0] * A[iz + 1*ixs - 1*iys - 1] + cx[0] * A[iz - 1*ixs + 1*iys + 1];
                    if constexpr(order >= 4)
                        bgx[iz] += -cx[1] * A[iz + 2*ixs - 2*iys - 2] + cx[1] * A[iz - 2*ixs + 2*iys + 2];
                    if constexpr(order >= 6)
                        bgx[iz] += -cx[2] * A[iz + 3*ixs - 3*iys - 3] + cx[2] * A[iz - 3*ixs + 3*iys + 3];
                    if constexpr(order >= 8)
                        bgx[iz] += -cx[3] * A[iz + 4*ixs - 4*iys - 4] + cx[3] * A[iz - 4*ixs + 4*iys + 4];
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

    double c1, c2 = 0.0;
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    double coeff0 = 0.0;
    if(order == 2) {c1=1.0;c2 = 0.0;}    // no optimzation for 2nd order
    for(int ax=0;ax < 13;ax++)
    {
        coeff0 += c1*LC2->plane_centers[ax] - c2*LC1->plane_centers[ax];
    }
    return scale*coeff0;
}


// Computes combined coefficients
template <typename RmgType>
void FiniteDiff::fd_combined_coeffs(int order, double hxgrid, int ax, RmgType * cm, RmgType *cp, double *kvec)
{
    double s1 = 2.0;
    RmgType t1, t2, t3, t4;
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

    double c1, c2=0.0;
    if(this->alt_laplacian && order > 2) c2 = cfac[0];
    c1 = scale*(1.0 + c2);
    t1 = c1*LC2->axis_lc[ax][3] - c2*LC1->axis_lc[ax][2];
    t2 = c1*LC2->axis_lc[ax][2] - c2*LC1->axis_lc[ax][1];
    t3 = c1*LC2->axis_lc[ax][1] - c2*LC1->axis_lc[ax][0];
    t4 = c1*LC2->axis_lc[ax][0];

    for(int i=0;i < order/2-1;i++)
    {
        t1 = c1*LC2->axis_lc[ax][order/2-i-1] - c2*LC1->axis_lc[ax][order/2-i-2];
        x1 = c1*LC2->axis_gc_x[ax][order/2-i-1] - c2*LC1->axis_gc_x[ax][order/2-i-2];
        y1 = c1*LC2->axis_gc_y[ax][order/2-i-1] - c2*LC1->axis_gc_y[ax][order/2-i-2];
        z1 = c1*LC2->axis_gc_z[ax][order/2-i-1] - c2*LC1->axis_gc_z[ax][order/2-i-2];
        cm[i] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
        cp[i] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    }
    t1 = c1*LC2->axis_lc[ax][0];
    x1 = c1*LC2->axis_gc_x[ax][0];
    y1 = c1*LC2->axis_gc_y[ax][0];
    z1 = c1*LC2->axis_gc_z[ax][0];
    cm[order/2-1] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[order/2-1] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

#if 0
    x1 = c1*LC2->axis_gc_x[ax][3] - c2*LC1->axis_gc_x[ax][2];
    y1 = c1*LC2->axis_gc_y[ax][3] - c2*LC1->axis_gc_y[ax][2];
    z1 = c1*LC2->axis_gc_z[ax][3] - c2*LC1->axis_gc_z[ax][2];
    cm[0] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[0] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC2->axis_gc_x[ax][2] - c2*LC1->axis_gc_x[ax][1];
    y1 = c1*LC2->axis_gc_y[ax][2] - c2*LC1->axis_gc_y[ax][1];
    z1 = c1*LC2->axis_gc_z[ax][2] - c2*LC1->axis_gc_z[ax][1];
    cm[1] = t2 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[1] = t2 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC2->axis_gc_x[ax][1] - c2*LC1->axis_gc_x[ax][0];
    y1 = c1*LC2->axis_gc_y[ax][1] - c2*LC1->axis_gc_y[ax][0];
    z1 = c1*LC2->axis_gc_z[ax][1] - c2*LC1->axis_gc_z[ax][0];
    cm[2] = t3 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[2] = t3 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC2->axis_gc_x[ax][0];
    y1 = c1*LC2->axis_gc_y[ax][0];
    z1 = c1*LC2->axis_gc_z[ax][0];
    cm[3] = t4 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
    cp[3] = t4 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);
#endif
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
    if(this->alt_laplacian && order > 2) c2 = cfac[0];
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
#if 0
    cx[0] = c1*LC2->axis_gc_x[axis][3] - c2*LC1->axis_gc_x[axis][2];
    cx[1] = c1*LC2->axis_gc_x[axis][2] - c2*LC1->axis_gc_x[axis][1];
    cx[2] = c1*LC2->axis_gc_x[axis][1] - c2*LC1->axis_gc_x[axis][0];
    cx[3] = c1*LC2->axis_gc_x[axis][0];

    cy[0] = c1*LC2->axis_gc_y[axis][3] - c2*LC1->axis_gc_y[axis][2];
    cy[1] = c1*LC2->axis_gc_y[axis][2] - c2*LC1->axis_gc_y[axis][1];
    cy[2] = c1*LC2->axis_gc_y[axis][1] - c2*LC1->axis_gc_y[axis][0];
    cy[3] = c1*LC2->axis_gc_y[axis][0];

    cz[0] = c1*LC2->axis_gc_z[axis][3] - c2*LC1->axis_gc_z[axis][2];
    cz[1] = c1*LC2->axis_gc_z[axis][2] - c2*LC1->axis_gc_z[axis][1];
    cz[2] = c1*LC2->axis_gc_z[axis][1] - c2*LC1->axis_gc_z[axis][0];
    cz[3] = c1*LC2->axis_gc_z[axis][0];
#endif
}


