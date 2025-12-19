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


// Handles case when all components are needed.
template <typename RmgType, int order>
void FiniteDiff::fd_gradient_general1 (RmgType * __restrict__ a, 
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
}

template void FiniteDiff::fd_gradient_general1<float, 6> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<double, 6> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<double>, 6> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<float>, 6> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general1<float, 8> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<double, 8> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<double>, 8> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<float>, 8> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general1<float, 10> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<double, 10> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<double>, 10> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<float>, 10> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);

template void FiniteDiff::fd_gradient_general1<float, 12> (float *, float *, float *, float *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<double, 12> (double *, double *, double *, double *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<double>, 12> (std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, std::complex<double>  *, double, int, int, int);
template void FiniteDiff::fd_gradient_general1<std::complex<float>, 12> (std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, std::complex<float>  *, double, int, int, int);


