/*
 *
 * Copyright 2024 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

/*
   This function is used in the TDDFT code to perform the reduction below on the GPU.

        for(st1 = 0; st1 < numst; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += psi[st1 * pbasis + idx] * xpsi[st1 * pbasis + idx];

*/

#if HIP_ENABLED


#include <complex>
#include <hip/hip_ext.h>
#include "Gpufuncs.h"

const int nTPB = 128; // Only 1 block but 128 threads per block


    template <typename T>
__global__ void vxcpsi(const std::complex<T> * __restrict__ psi, std::complex<T> * __restrict__ vpsi, const T * __restrict__ vxc_psi, const int pbasis, const int numst){

  //for(size_t pidx = 0;pidx < pbasis;pidx += nTPB)
    size_t pidx = blockIdx.x*nTPB;
    int idx = pidx + threadIdx.x;
    if(idx < pbasis)
    {
        for(size_t ns = 0;ns < numst;ns++)
        {
            vpsi[idx] += psi[idx] * std::complex<T>(vxc_psi[3*pbasis+idx], 0.0);
            vpsi[idx] += psi[idx+pbasis] * std::complex<T>(vxc_psi[pbasis+idx], vxc_psi[2*pbasis+idx]);
            vpsi[idx + pbasis] += - psi[idx + pbasis] * std::complex<T>(vxc_psi[3*pbasis+idx], 0.0);
            vpsi[idx + pbasis] += psi[idx] * std::complex<T>(vxc_psi[pbasis+idx], -vxc_psi[2*pbasis+idx]);

        }
    }
}

template void GpuVxc_x_psi_noncoll<double>(std::complex<double> *psi, std::complex<double> *xpsi, double *vxc, int pbasis, int num_states);
template void GpuVxc_x_psi_noncoll<float>(std::complex<float> *psi, std::complex<float> *xpsi, float *vxc, int pbasis, int num_states);
    template <typename TypeV>
void GpuVxc_x_psi_noncoll(std::complex<TypeV> *psi, std::complex<TypeV> *xpsi, TypeV *vxc, int pbasis, int num_states)
{
    int nblocks = pbasis / nTPB + 1;
    //  vxc_noncoll * psi (incuding 2 spinors)
    vxcpsi<<<nblocks, nTPB>>>(psi, xpsi, vxc, pbasis, num_states);
}


#endif
