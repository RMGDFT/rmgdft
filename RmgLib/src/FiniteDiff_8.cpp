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

template void FiniteDiff::app8_gradient_general<float> (float *, float *, float *, float *, int, int, int);
template void FiniteDiff::app8_gradient_general<double> (double *, double *, double *, double *, int, int, int);
template void FiniteDiff::app8_gradient_general<std::complex<double> > (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, int, int, int);
template void FiniteDiff::app8_gradient_general<std::complex<float> > (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, int, int, int);

template double FiniteDiff::app8_combined<float>(float *, float *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<double>(double *, double *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<std::complex <float> >(std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, double *kvec, bool use_gpu);
template double FiniteDiff::app8_combined<std::complex <double> >(std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, double *kvec, bool use_gpu);

template void FiniteDiff::app8_gradient_coeffs<float>(int , int , float *, float *, float *);
template void FiniteDiff::app8_gradient_coeffs<double>(int , int , double *, double *, double *);
template void FiniteDiff::app8_gradient_coeffs<std::complex<double>>(int , int , std::complex<double> *, std::complex<double> *, std::complex<double> *);

template void FiniteDiff::app8_combined_coeffs<float>(int , int , float *, float *, double *);
template void FiniteDiff::app8_combined_coeffs<double>(int , int , double *, double *, double *);
template void FiniteDiff::app8_combined_coeffs<std::complex<float>>(int , int , std::complex<float> *, std::complex<float> *, double *);
template void FiniteDiff::app8_combined_coeffs<std::complex<double>>(int , int , std::complex<double> *, std::complex<double> *, double *);


template <typename RmgType>
double FiniteDiff::app8_del2(RmgType * __restrict__ a, RmgType * __restrict__ b, int dimx, int dimy, int dimz,
               double gridhx, double gridhy, double gridhz)
{
    double kvec[3] = {0.0,0.0,0.0};
    return FiniteDiff::app8_combined(a, b, dimx, dimy, dimz,
        gridhx, gridhy, gridhz, kvec, false);

}  /* end app8_del2 */


template <typename RmgType>
void FiniteDiff::app_gradient_eighth (RmgType * __restrict__ rptr, RmgType * __restrict__ wxr, RmgType * __restrict__ wyr, RmgType * __restrict__ wzr, int dimx, int dimy, int dimz,
        double gridhx, double gridhy, double gridhz)
{
    FiniteDiff::app8_gradient_general (rptr, wxr, wyr, wzr, dimx, dimy, dimz);
}



template <typename RmgType>
double FiniteDiff::app8_combined(RmgType * __restrict__ a, RmgType * __restrict__ b, 
		int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
		double *kvec, bool use_gpu)
{
    int ibrav = L->get_ibrav_type();
    RmgType cm[4], cp[4];
    RmgType cpx[4], cmx[4], cpy[4], cmy[4], cpz[4], cmz[4];
    int ixs = (dimy + 8) * (dimz + 8);
    int iys = (dimz + 8);


    // NULL b means we just want the diagonal component.
    double th2 = app8_coeff0();
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
    app8_combined_coeffs(8, 0, cmx, cpx, kvec);
    app8_combined_coeffs(8, 1, cmy, cpy, kvec);
    app8_combined_coeffs(8, 2, cmz, cpz, kvec);

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
                    cpz[0] * A[iz + 1] + cmz[0] * A[iz - 1] +
                    cpz[1] * A[iz + 2] + cmz[1] * A[iz - 2] +
                    cpz[2] * A[iz + 3] + cmz[2] * A[iz - 3] +
                    cpz[3] * A[iz + 4] + cmz[3] * A[iz - 4];
            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    cpy[0] * A[iz + iys] + cmy[0] * A[iz - iys] +
                    cpy[1] * A[iz + 2*iys] + cmy[1] * A[iz - 2*iys] +
                    cpy[2] * A[iz + 3*iys] + cmy[2] * A[iz - 3*iys] +
                    cpy[3] * A[iz + 4*iys] + cmy[3] * A[iz - 4*iys];
            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                B[iz] +=
                    cpx[0] * A[iz + ixs] + cmx[0] * A[iz - ixs] +
                    cpx[1] * A[iz + 2*ixs] + cmx[1] * A[iz - 2*ixs] +
                    cpx[2] * A[iz + 3*ixs] + cmx[2] * A[iz - 3*ixs] +
                    cpx[3] * A[iz + 4*ixs] + cmx[3] * A[iz - 4*ixs];
            }                   /* end for */
        }
    }

    /* Quick return for orthogonal axis cases */
    if(ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == CUBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
        return (double)std::real(th2);

    // Add additional axes as required
    if(LC->include_axis[3])
    {
        app8_combined_coeffs(8, 3, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + ixs + iys] + cm[0] * A[iz - ixs - iys] +
                        cp[1] * A[iz + 2*ixs + 2*iys] + cm[1] * A[iz - 2*ixs - 2*iys] +
                        cp[2] * A[iz + 3*ixs + 3*iys] + cm[2] * A[iz - 3*ixs - 3*iys] +
                        cp[3] * A[iz + 4*ixs + 4*iys] + cm[3] * A[iz - 4*ixs - 4*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[4])
    {
        app8_combined_coeffs(8, 4, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + ixs + 1] + cm[0] * A[iz - ixs - 1] +
                        cp[1] * A[iz + 2*ixs + 2] + cm[1] * A[iz - 2*ixs - 2] +
                        cp[2] * A[iz + 3*ixs + 3] + cm[2] * A[iz - 3*ixs - 3] +
                        cp[3] * A[iz + 4*ixs + 4] + cm[3] * A[iz - 4*ixs - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[5])
    {
        app8_combined_coeffs(8, 5, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + iys + 1] + cm[0] * A[iz - iys - 1] +
                        cp[1] * A[iz + 2*iys + 2] + cm[1] * A[iz - 2*iys - 2] +
                        cp[2] * A[iz + 3*iys + 3] + cm[2] * A[iz - 3*iys - 3] +
                        cp[3] * A[iz + 4*iys + 4] + cm[3] * A[iz - 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[6])
    {
        app8_combined_coeffs(8, 6, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz - ixs + iys] + cm[0] * A[iz + ixs - iys] +
                        cp[1] * A[iz - 2*ixs + 2*iys] + cm[1] * A[iz + 2*ixs - 2*iys] +
                        cp[2] * A[iz - 3*ixs + 3*iys] + cm[2] * A[iz + 3*ixs - 3*iys] +
                        cp[3] * A[iz - 4*ixs + 4*iys] + cm[3] * A[iz + 4*ixs - 4*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[7])
    {
        app8_combined_coeffs(8, 7, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz - ixs + 1] + cm[0] * A[iz + ixs - 1] +
                        cp[1] * A[iz - 2*ixs + 2] + cm[1] * A[iz + 2*ixs - 2] +
                        cp[2] * A[iz - 3*ixs + 3] + cm[2] * A[iz + 3*ixs - 3] +
                        cp[3] * A[iz - 4*ixs + 4] + cm[3] * A[iz + 4*ixs - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[8])
    {
        app8_combined_coeffs(8, 8, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz - iys + 1] + cm[0] * A[iz + iys - 1] +
                        cp[1] * A[iz - 2*iys + 2] + cm[1] * A[iz + 2*iys - 2] +
                        cp[2] * A[iz - 3*iys + 3] + cm[2] * A[iz + 3*iys - 3] +
                        cp[3] * A[iz - 4*iys + 4] + cm[3] * A[iz + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[9])
    {
        app8_combined_coeffs(8, 9, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + 1*ixs + 1*iys + 1] + cm[0] * A[iz - 1*ixs - 1*iys - 1] +
                        cp[1] * A[iz + 2*ixs + 2*iys + 2] + cm[1] * A[iz - 2*ixs - 2*iys - 2] +
                        cp[2] * A[iz + 3*ixs + 3*iys + 3] + cm[2] * A[iz - 3*ixs - 3*iys - 3] +
                        cp[3] * A[iz + 4*ixs + 4*iys + 4] + cm[3] * A[iz - 4*ixs - 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[10])
    {
        app8_combined_coeffs(8, 10, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz - 1*ixs - 1*iys + 1] + cm[0] * A[iz + 1*ixs + 1*iys - 1] +
                        cp[1] * A[iz - 2*ixs - 2*iys + 2] + cm[1] * A[iz + 2*ixs + 2*iys - 2] +
                        cp[2] * A[iz - 3*ixs - 3*iys + 3] + cm[2] * A[iz + 3*ixs + 3*iys - 3] +
                        cp[3] * A[iz - 4*ixs - 4*iys + 4] + cm[3] * A[iz + 4*ixs + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[11])
    {
        app8_combined_coeffs(8, 11, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + 1*ixs - 1*iys + 1] + cm[0] * A[iz - 1*ixs + 1*iys - 1] +
                        cp[1] * A[iz + 2*ixs - 2*iys + 2] + cm[1] * A[iz - 2*ixs + 2*iys - 2] +
                        cp[2] * A[iz + 3*ixs - 3*iys + 3] + cm[2] * A[iz - 3*ixs + 3*iys - 3] +
                        cp[3] * A[iz + 4*ixs - 4*iys + 4] + cm[3] * A[iz - 4*ixs + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[12])
    {
        app8_combined_coeffs(8, 12, cm, cp, kvec);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *B = &b[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    B[iz] +=
                        cp[0] * A[iz + 1*ixs - 1*iys - 1] + cm[0] * A[iz - 1*ixs + 1*iys + 1] +
                        cp[1] * A[iz + 2*ixs - 2*iys - 2] + cm[1] * A[iz - 2*ixs + 2*iys + 2] +
                        cp[2] * A[iz + 3*ixs - 3*iys - 3] + cm[2] * A[iz - 3*ixs + 3*iys + 3] +
                        cp[3] * A[iz + 4*ixs - 4*iys - 4] + cm[3] * A[iz - 4*ixs + 4*iys + 4];
                }                   /* end for */
            }
        }
    }

    /* Return the diagonal component of the operator */
    return (double)std::real(th2);

} /* end app8_combined */

// Gets the central coefficient
double FiniteDiff::app8_coeff0(void)
{
    double c1, c2 = 0.0;
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    double coeff0 = 0.0;
    for(int ax=0;ax < 13;ax++)
    {
        coeff0 += c1*LC->plane_centers[ax] - c2*LC_6->plane_centers[ax];
    }
    return coeff0;
}


// Computes combined coefficients
template <typename RmgType>
void FiniteDiff::app8_combined_coeffs(int order, int ax, RmgType * cm, RmgType *cp, double *kvec)
{
    double s1 = 2.0;
    RmgType t1, t2, t3, t4;


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

    double c1, c2=0.0;
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    t1 = c1*LC->axis_lc[ax][3] - c2*LC_6->axis_lc[ax][2];
    t2 = c1*LC->axis_lc[ax][2] - c2*LC_6->axis_lc[ax][1];
    t3 = c1*LC->axis_lc[ax][1] - c2*LC_6->axis_lc[ax][0];
    t4 = c1*LC->axis_lc[ax][0];

    RmgType x1, y1, z1;
    x1 = c1*LC->axis_gc_x[ax][3] - c2*LC_6->axis_gc_x[ax][2];
    y1 = c1*LC->axis_gc_y[ax][3] - c2*LC_6->axis_gc_y[ax][2];
    z1 = c1*LC->axis_gc_z[ax][3] - c2*LC_6->axis_gc_z[ax][2];
    cm[0] = t1 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][2] - c2*LC_6->axis_gc_x[ax][1];
    y1 = c1*LC->axis_gc_y[ax][2] - c2*LC_6->axis_gc_y[ax][1];
    z1 = c1*LC->axis_gc_z[ax][2] - c2*LC_6->axis_gc_z[ax][1];
    cm[1] = t2 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][1] - c2*LC_6->axis_gc_x[ax][0];
    y1 = c1*LC->axis_gc_y[ax][1] - c2*LC_6->axis_gc_y[ax][0];
    z1 = c1*LC->axis_gc_z[ax][1] - c2*LC_6->axis_gc_z[ax][0];
    cm[2] = t3 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][0];
    y1 = c1*LC->axis_gc_y[ax][0];
    z1 = c1*LC->axis_gc_z[ax][0];
    cm[3] = t4 + s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][3] - c2*LC_6->axis_gc_x[ax][2];
    y1 = c1*LC->axis_gc_y[ax][3] - c2*LC_6->axis_gc_y[ax][2];
    z1 = c1*LC->axis_gc_z[ax][3] - c2*LC_6->axis_gc_z[ax][2];
    cp[0] = t1 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][2] - c2*LC_6->axis_gc_x[ax][1];
    y1 = c1*LC->axis_gc_y[ax][2] - c2*LC_6->axis_gc_y[ax][1];
    z1 = c1*LC->axis_gc_z[ax][2] - c2*LC_6->axis_gc_z[ax][1];
    cp[1] = t2 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][1] - c2*LC_6->axis_gc_x[ax][0];
    y1 = c1*LC->axis_gc_y[ax][1] - c2*LC_6->axis_gc_y[ax][0];
    z1 = c1*LC->axis_gc_z[ax][1] - c2*LC_6->axis_gc_z[ax][0];
    cp[2] = t3 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

    x1 = c1*LC->axis_gc_x[ax][0];
    y1 = c1*LC->axis_gc_y[ax][0];
    z1 = c1*LC->axis_gc_z[ax][0];
    cp[3] = t4 - s1 * I_t * (kvec[0]*x1 + kvec[1]*y1 + kvec[2]*z1);

}

template <typename RmgType>
void FiniteDiff::app8_gradient_coeffs(int order, int axis , RmgType *cx, RmgType *cy, RmgType *cz)
{
    double c1, c2=0.0;
    if(this->alt_laplacian) c2 = cfac[0];
    c1 = 1.0 + c2;
    cx[0] = c1*LC->axis_gc_x[axis][3] - c2*LC_6->axis_gc_x[axis][2];
    cx[1] = c1*LC->axis_gc_x[axis][2] - c2*LC_6->axis_gc_x[axis][1];
    cx[2] = c1*LC->axis_gc_x[axis][1] - c2*LC_6->axis_gc_x[axis][0];
    cx[3] = c1*LC->axis_gc_x[axis][0];

    cy[0] = c1*LC->axis_gc_y[axis][3] - c2*LC_6->axis_gc_y[axis][2];
    cy[1] = c1*LC->axis_gc_y[axis][2] - c2*LC_6->axis_gc_y[axis][1];
    cy[2] = c1*LC->axis_gc_y[axis][1] - c2*LC_6->axis_gc_y[axis][0];
    cy[3] = c1*LC->axis_gc_y[axis][0];

    cz[0] = c1*LC->axis_gc_z[axis][3] - c2*LC_6->axis_gc_z[axis][2];
    cz[1] = c1*LC->axis_gc_z[axis][2] - c2*LC_6->axis_gc_z[axis][1];
    cz[2] = c1*LC->axis_gc_z[axis][1] - c2*LC_6->axis_gc_z[axis][0];
    cz[3] = c1*LC->axis_gc_z[axis][0];
}



template <typename RmgType>
void FiniteDiff::app8_gradient_general (RmgType * __restrict__ a, 
                                RmgType * __restrict__ gx, 
                                RmgType * __restrict__ gy, 
                                RmgType * __restrict__ gz,
                                int dimx, int dimy, int dimz)
{
    int ibrav = L->get_ibrav_type();
    RmgType cxx[4], cxy[4], cxz[4];
    RmgType cyx[4], cyy[4], cyz[4];
    RmgType czx[4], czy[4], czz[4];
    RmgType cx[4], cy[4], cz[4];
    int ixs = (dimy + 8) * (dimz + 8);
    int iys = (dimz + 8);


    // Get coeffs for x,y,z axes which are used by all lattice types
    app8_gradient_coeffs(8, 0, cxx, cxy, cxz);
    app8_gradient_coeffs(8, 1, cyx, cyy, cyz);
    app8_gradient_coeffs(8, 2, czx, czy, czz);
    for (int ix = 4; ix < dimx + 4; ix++)
    {
        for (int iy = 4; iy < dimy + 4; iy++)
        {
            RmgType *A = &a[iy*iys + ix*ixs];
            RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
            RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
            RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
            // z-direction is orthogonal to xy-plane and only requires increments/decrements along z
            // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                bgx[iz] =
                    -czx[0] * A[iz + 1] + czx[0] * A[iz - 1] +
                    -czx[1] * A[iz + 2] + czx[1] * A[iz - 2] +
                    -czx[2] * A[iz + 3] + czx[2] * A[iz - 3] +
                    -czx[3] * A[iz + 4] + czx[3] * A[iz - 4];
                bgy[iz] =
                    -czy[0] * A[iz + 1] + czy[0] * A[iz - 1] +
                    -czy[1] * A[iz + 2] + czy[1] * A[iz - 2] +
                    -czy[2] * A[iz + 3] + czy[2] * A[iz - 3] +
                    -czy[3] * A[iz + 4] + czy[3] * A[iz - 4];
                bgz[iz] =
                    -czz[0] * A[iz + 1] + czz[0] * A[iz - 1] +
                    -czz[1] * A[iz + 2] + czz[1] * A[iz - 2] +
                    -czz[2] * A[iz + 3] + czz[2] * A[iz - 3] +
                    -czz[3] * A[iz + 4] + czz[3] * A[iz - 4];

            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                bgx[iz] +=
                    -cyx[0] * A[iz + iys] + cyx[0] * A[iz - iys] +
                    -cyx[1] * A[iz + 2*iys] + cyx[1] * A[iz - 2*iys] +
                    -cyx[2] * A[iz + 3*iys] + cyx[2] * A[iz - 3*iys] +
                    -cyx[3] * A[iz + 4*iys] + cyx[3] * A[iz - 4*iys];
                bgy[iz] +=
                    -cyy[0] * A[iz + iys] + cyy[0] * A[iz - iys] +
                    -cyy[1] * A[iz + 2*iys] + cyy[1] * A[iz - 2*iys] +
                    -cyy[2] * A[iz + 3*iys] + cyy[2] * A[iz - 3*iys] +
                    -cyy[3] * A[iz + 4*iys] + cyy[3] * A[iz - 4*iys];
                bgz[iz] +=
                    -cyz[0] * A[iz + iys] + cyz[0] * A[iz - iys] +
                    -cyz[1] * A[iz + 2*iys] + cyz[1] * A[iz - 2*iys] +
                    -cyz[2] * A[iz + 3*iys] + cyz[2] * A[iz - 3*iys] +
                    -cyz[3] * A[iz + 4*iys] + cyz[3] * A[iz - 4*iys];
            }
            for (int iz = 4; iz < dimz + 4; iz++)
            {
                bgx[iz] +=
                    -cxx[0] * A[iz + ixs] + cxx[0] * A[iz - ixs] +
                    -cxx[1] * A[iz + 2*ixs] + cxx[1] * A[iz - 2*ixs] +
                    -cxx[2] * A[iz + 3*ixs] + cxx[2] * A[iz - 3*ixs] +
                    -cxx[3] * A[iz + 4*ixs] + cxx[3] * A[iz - 4*ixs];
                bgy[iz] +=
                    -cxy[0] * A[iz + ixs] + cxy[0] * A[iz - ixs] +
                    -cxy[1] * A[iz + 2*ixs] + cxy[1] * A[iz - 2*ixs] +
                    -cxy[2] * A[iz + 3*ixs] + cxy[2] * A[iz - 3*ixs] +
                    -cxy[3] * A[iz + 4*ixs] + cxy[3] * A[iz - 4*ixs];
                bgz[iz] +=
                    -cxz[0] * A[iz + ixs] + cxz[0] * A[iz - ixs] +
                    -cxz[1] * A[iz + 2*ixs] + cxz[1] * A[iz - 2*ixs] +
                    -cxz[2] * A[iz + 3*ixs] + cxz[2] * A[iz - 3*ixs] +
                    -cxz[3] * A[iz + 4*ixs] + cxz[3] * A[iz - 4*ixs];

            }                   /* end for */
        }
    }

    /* Quick return for orthogonal axis cases */
    if(ibrav == ORTHORHOMBIC_PRIMITIVE || ibrav == CUBIC_PRIMITIVE || ibrav == TETRAGONAL_PRIMITIVE)
        return;

    // Add additional axes as required
    if(LC->include_axis[3])
    {
        app8_gradient_coeffs(8, 3, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + ixs + iys] + cz[0] * A[iz - ixs - iys] +
                        -cz[1] * A[iz + 2*ixs + 2*iys] + cz[1] * A[iz - 2*ixs - 2*iys] +
                        -cz[2] * A[iz + 3*ixs + 3*iys] + cz[2] * A[iz - 3*ixs - 3*iys] +
                        -cz[3] * A[iz + 4*ixs + 4*iys] + cz[3] * A[iz - 4*ixs - 4*iys];
                    bgy[iz] +=
                        -cy[0] * A[iz + ixs + iys] + cy[0] * A[iz - ixs - iys] +
                        -cy[1] * A[iz + 2*ixs + 2*iys] + cy[1] * A[iz - 2*ixs - 2*iys] +
                        -cy[2] * A[iz + 3*ixs + 3*iys] + cy[2] * A[iz - 3*ixs - 3*iys] +
                        -cy[3] * A[iz + 4*ixs + 4*iys] + cy[3] * A[iz - 4*ixs - 4*iys];
                    bgx[iz] +=
                        -cx[0] * A[iz + ixs + iys] + cx[0] * A[iz - ixs - iys] +
                        -cx[1] * A[iz + 2*ixs + 2*iys] + cx[1] * A[iz - 2*ixs - 2*iys] +
                        -cx[2] * A[iz + 3*ixs + 3*iys] + cx[2] * A[iz - 3*ixs - 3*iys] +
                        -cx[3] * A[iz + 4*ixs + 4*iys] + cx[3] * A[iz - 4*ixs - 4*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[4])
    {
        app8_gradient_coeffs(8, 4, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + ixs + 1] + cz[0] * A[iz - ixs - 1] +
                        -cz[1] * A[iz + 2*ixs + 2] + cz[1] * A[iz - 2*ixs - 2] +
                        -cz[2] * A[iz + 3*ixs + 3] + cz[2] * A[iz - 3*ixs - 3] +
                        -cz[3] * A[iz + 4*ixs + 4] + cz[3] * A[iz - 4*ixs - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz + ixs + 1] + cy[0] * A[iz - ixs - 1] +
                        -cy[1] * A[iz + 2*ixs + 2] + cy[1] * A[iz - 2*ixs - 2] +
                        -cy[2] * A[iz + 3*ixs + 3] + cy[2] * A[iz - 3*ixs - 3] +
                        -cy[3] * A[iz + 4*ixs + 4] + cy[3] * A[iz - 4*ixs - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz + ixs + 1] + cx[0] * A[iz - ixs - 1] +
                        -cx[1] * A[iz + 2*ixs + 2] + cx[1] * A[iz - 2*ixs - 2] +
                        -cx[2] * A[iz + 3*ixs + 3] + cx[2] * A[iz - 3*ixs - 3] +
                        -cx[3] * A[iz + 4*ixs + 4] + cx[3] * A[iz - 4*ixs - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[5])
    {
        app8_gradient_coeffs(8, 5, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + iys + 1] + cz[0] * A[iz - iys - 1] +
                        -cz[1] * A[iz + 2*iys + 2] + cz[1] * A[iz - 2*iys - 2] +
                        -cz[2] * A[iz + 3*iys + 3] + cz[2] * A[iz - 3*iys - 3] +
                        -cz[3] * A[iz + 4*iys + 4] + cz[3] * A[iz - 4*iys - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz + iys + 1] + cy[0] * A[iz - iys - 1] +
                        -cy[1] * A[iz + 2*iys + 2] + cy[1] * A[iz - 2*iys - 2] +
                        -cy[2] * A[iz + 3*iys + 3] + cy[2] * A[iz - 3*iys - 3] +
                        -cy[3] * A[iz + 4*iys + 4] + cy[3] * A[iz - 4*iys - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz + iys + 1] + cx[0] * A[iz - iys - 1] +
                        -cx[1] * A[iz + 2*iys + 2] + cx[1] * A[iz - 2*iys - 2] +
                        -cx[2] * A[iz + 3*iys + 3] + cx[2] * A[iz - 3*iys - 3] +
                        -cx[3] * A[iz + 4*iys + 4] + cx[3] * A[iz - 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[6])
    {
        app8_gradient_coeffs(8, 6, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz - ixs + iys] + cz[0] * A[iz + ixs - iys] +
                        -cz[1] * A[iz - 2*ixs + 2*iys] + cz[1] * A[iz + 2*ixs - 2*iys] +
                        -cz[2] * A[iz - 3*ixs + 3*iys] + cz[2] * A[iz + 3*ixs - 3*iys] +
                        -cz[3] * A[iz - 4*ixs + 4*iys] + cz[3] * A[iz + 4*ixs - 4*iys];
                    bgy[iz] +=
                        -cy[0] * A[iz - ixs + iys] + cy[0] * A[iz + ixs - iys] +
                        -cy[1] * A[iz - 2*ixs + 2*iys] + cy[1] * A[iz + 2*ixs - 2*iys] +
                        -cy[2] * A[iz - 3*ixs + 3*iys] + cy[2] * A[iz + 3*ixs - 3*iys] +
                        -cy[3] * A[iz - 4*ixs + 4*iys] + cy[3] * A[iz + 4*ixs - 4*iys];
                    bgx[iz] +=
                        -cx[0] * A[iz - ixs + iys] + cx[0] * A[iz + ixs - iys] +
                        -cx[1] * A[iz - 2*ixs + 2*iys] + cx[1] * A[iz + 2*ixs - 2*iys] +
                        -cx[2] * A[iz - 3*ixs + 3*iys] + cx[2] * A[iz + 3*ixs - 3*iys] +
                        -cx[3] * A[iz - 4*ixs + 4*iys] + cx[3] * A[iz + 4*ixs - 4*iys];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[7])
    {
        app8_gradient_coeffs(8, 7, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz - ixs + 1] + cz[0] * A[iz + ixs - 1] +
                        -cz[1] * A[iz - 2*ixs + 2] + cz[1] * A[iz + 2*ixs - 2] +
                        -cz[2] * A[iz - 3*ixs + 3] + cz[2] * A[iz + 3*ixs - 3] +
                        -cz[3] * A[iz - 4*ixs + 4] + cz[3] * A[iz + 4*ixs - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz - ixs + 1] + cy[0] * A[iz + ixs - 1] +
                        -cy[1] * A[iz - 2*ixs + 2] + cy[1] * A[iz + 2*ixs - 2] +
                        -cy[2] * A[iz - 3*ixs + 3] + cy[2] * A[iz + 3*ixs - 3] +
                        -cy[3] * A[iz - 4*ixs + 4] + cy[3] * A[iz + 4*ixs - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz - ixs + 1] + cx[0] * A[iz + ixs - 1] +
                        -cx[1] * A[iz - 2*ixs + 2] + cx[1] * A[iz + 2*ixs - 2] +
                        -cx[2] * A[iz - 3*ixs + 3] + cx[2] * A[iz + 3*ixs - 3] +
                        -cx[3] * A[iz - 4*ixs + 4] + cx[3] * A[iz + 4*ixs - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[8])
    {
        app8_gradient_coeffs(8, 8, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz - iys + 1] + cz[0] * A[iz + iys - 1] +
                        -cz[1] * A[iz - 2*iys + 2] + cz[1] * A[iz + 2*iys - 2] +
                        -cz[2] * A[iz - 3*iys + 3] + cz[2] * A[iz + 3*iys - 3] +
                        -cz[3] * A[iz - 4*iys + 4] + cz[3] * A[iz + 4*iys - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz - iys + 1] + cy[0] * A[iz + iys - 1] +
                        -cy[1] * A[iz - 2*iys + 2] + cy[1] * A[iz + 2*iys - 2] +
                        -cy[2] * A[iz - 3*iys + 3] + cy[2] * A[iz + 3*iys - 3] +
                        -cy[3] * A[iz - 4*iys + 4] + cy[3] * A[iz + 4*iys - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz - iys + 1] + cx[0] * A[iz + iys - 1] +
                        -cx[1] * A[iz - 2*iys + 2] + cx[1] * A[iz + 2*iys - 2] +
                        -cx[2] * A[iz - 3*iys + 3] + cx[2] * A[iz + 3*iys - 3] +
                        -cx[3] * A[iz - 4*iys + 4] + cx[3] * A[iz + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[9])
    {
        app8_gradient_coeffs(8, 9, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + 1*ixs + 1*iys + 1] + cz[0] * A[iz - 1*ixs - 1*iys - 1] +
                        -cz[1] * A[iz + 2*ixs + 2*iys + 2] + cz[1] * A[iz - 2*ixs - 2*iys - 2] +
                        -cz[2] * A[iz + 3*ixs + 3*iys + 3] + cz[2] * A[iz - 3*ixs - 3*iys - 3] +
                        -cz[3] * A[iz + 4*ixs + 4*iys + 4] + cz[3] * A[iz - 4*ixs - 4*iys - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz + 1*ixs + 1*iys + 1] + cy[0] * A[iz - 1*ixs - 1*iys - 1] +
                        -cy[1] * A[iz + 2*ixs + 2*iys + 2] + cy[1] * A[iz - 2*ixs - 2*iys - 2] +
                        -cy[2] * A[iz + 3*ixs + 3*iys + 3] + cy[2] * A[iz - 3*ixs - 3*iys - 3] +
                        -cy[3] * A[iz + 4*ixs + 4*iys + 4] + cy[3] * A[iz - 4*ixs - 4*iys - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz + 1*ixs + 1*iys + 1] + cx[0] * A[iz - 1*ixs - 1*iys - 1] +
                        -cx[1] * A[iz + 2*ixs + 2*iys + 2] + cx[1] * A[iz - 2*ixs - 2*iys - 2] +
                        -cx[2] * A[iz + 3*ixs + 3*iys + 3] + cx[2] * A[iz - 3*ixs - 3*iys - 3] +
                        -cx[3] * A[iz + 4*ixs + 4*iys + 4] + cx[3] * A[iz - 4*ixs - 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[10])
    {
        app8_gradient_coeffs(8, 10, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz - 1*ixs - 1*iys + 1] + cz[0] * A[iz + 1*ixs + 1*iys - 1] +
                        -cz[1] * A[iz - 2*ixs - 2*iys + 2] + cz[1] * A[iz + 2*ixs + 2*iys - 2] +
                        -cz[2] * A[iz - 3*ixs - 3*iys + 3] + cz[2] * A[iz + 3*ixs + 3*iys - 3] +
                        -cz[3] * A[iz - 4*ixs - 4*iys + 4] + cz[3] * A[iz + 4*ixs + 4*iys - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz - 1*ixs - 1*iys + 1] + cy[0] * A[iz + 1*ixs + 1*iys - 1] +
                        -cy[1] * A[iz - 2*ixs - 2*iys + 2] + cy[1] * A[iz + 2*ixs + 2*iys - 2] +
                        -cy[2] * A[iz - 3*ixs - 3*iys + 3] + cy[2] * A[iz + 3*ixs + 3*iys - 3] +
                        -cy[3] * A[iz - 4*ixs - 4*iys + 4] + cy[3] * A[iz + 4*ixs + 4*iys - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz - 1*ixs - 1*iys + 1] + cx[0] * A[iz + 1*ixs + 1*iys - 1] +
                        -cx[1] * A[iz - 2*ixs - 2*iys + 2] + cx[1] * A[iz + 2*ixs + 2*iys - 2] +
                        -cx[2] * A[iz - 3*ixs - 3*iys + 3] + cx[2] * A[iz + 3*ixs + 3*iys - 3] +
                        -cx[3] * A[iz - 4*ixs - 4*iys + 4] + cx[3] * A[iz + 4*ixs + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[11])
    {
        app8_gradient_coeffs(8, 11, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + 1*ixs - 1*iys + 1] + cz[0] * A[iz - 1*ixs + 1*iys - 1] +
                        -cz[1] * A[iz + 2*ixs - 2*iys + 2] + cz[1] * A[iz - 2*ixs + 2*iys - 2] +
                        -cz[2] * A[iz + 3*ixs - 3*iys + 3] + cz[2] * A[iz - 3*ixs + 3*iys - 3] +
                        -cz[3] * A[iz + 4*ixs - 4*iys + 4] + cz[3] * A[iz - 4*ixs + 4*iys - 4];
                    bgy[iz] +=
                        -cy[0] * A[iz + 1*ixs - 1*iys + 1] + cy[0] * A[iz - 1*ixs + 1*iys - 1] +
                        -cy[1] * A[iz + 2*ixs - 2*iys + 2] + cy[1] * A[iz - 2*ixs + 2*iys - 2] +
                        -cy[2] * A[iz + 3*ixs - 3*iys + 3] + cy[2] * A[iz - 3*ixs + 3*iys - 3] +
                        -cy[3] * A[iz + 4*ixs - 4*iys + 4] + cy[3] * A[iz - 4*ixs + 4*iys - 4];
                    bgx[iz] +=
                        -cx[0] * A[iz + 1*ixs - 1*iys + 1] + cx[0] * A[iz - 1*ixs + 1*iys - 1] +
                        -cx[1] * A[iz + 2*ixs - 2*iys + 2] + cx[1] * A[iz - 2*ixs + 2*iys - 2] +
                        -cx[2] * A[iz + 3*ixs - 3*iys + 3] + cx[2] * A[iz - 3*ixs + 3*iys - 3] +
                        -cx[3] * A[iz + 4*ixs - 4*iys + 4] + cx[3] * A[iz - 4*ixs + 4*iys - 4];
                }                   /* end for */
            }
        }
    }

    if(LC->include_axis[12])
    {
        app8_gradient_coeffs(8, 12, cx, cy, cz);
        for (int ix = 4; ix < dimx + 4; ix++)
        {
            for (int iy = 4; iy < dimy + 4; iy++)
            {
                RmgType *A = &a[iy*iys + ix*ixs];
                RmgType *bgx = &gx[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgy = &gy[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                RmgType *bgz = &gz[(iy - 4)*dimz + (ix - 4)*dimy*dimz - 4];
                for (int iz = 4; iz < dimz + 4; iz++)
                {
                    bgz[iz] +=
                        -cz[0] * A[iz + 1*ixs - 1*iys - 1] + cz[0] * A[iz - 1*ixs + 1*iys + 1] +
                        -cz[1] * A[iz + 2*ixs - 2*iys - 2] + cz[1] * A[iz - 2*ixs + 2*iys + 2] +
                        -cz[2] * A[iz + 3*ixs - 3*iys - 3] + cz[2] * A[iz - 3*ixs + 3*iys + 3] +
                        -cz[3] * A[iz + 4*ixs - 4*iys - 4] + cz[3] * A[iz - 4*ixs + 4*iys + 4];
                    bgy[iz] +=
                        -cy[0] * A[iz + 1*ixs - 1*iys - 1] + cy[0] * A[iz - 1*ixs + 1*iys + 1] +
                        -cy[1] * A[iz + 2*ixs - 2*iys - 2] + cy[1] * A[iz - 2*ixs + 2*iys + 2] +
                        -cy[2] * A[iz + 3*ixs - 3*iys - 3] + cy[2] * A[iz - 3*ixs + 3*iys + 3] +
                        -cy[3] * A[iz + 4*ixs - 4*iys - 4] + cy[3] * A[iz - 4*ixs + 4*iys + 4];
                    bgx[iz] +=
                        -cx[0] * A[iz + 1*ixs - 1*iys - 1] + cx[0] * A[iz - 1*ixs + 1*iys + 1] +
                        -cx[1] * A[iz + 2*ixs - 2*iys - 2] + cx[1] * A[iz - 2*ixs + 2*iys + 2] +
                        -cx[2] * A[iz + 3*ixs - 3*iys - 3] + cx[2] * A[iz - 3*ixs + 3*iys + 3] +
                        -cx[3] * A[iz + 4*ixs - 4*iys - 4] + cx[3] * A[iz - 4*ixs + 4*iys + 4];
                }                   /* end for */
            }
        }
    }
} /* end app8_gradient_general */

