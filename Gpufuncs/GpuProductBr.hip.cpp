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



// produce the sum or mean of each array
template <typename T>
__global__ void breduce(const T * __restrict__ idata1, const T * __restrict__ idata2, double * __restrict__ odata, const int numst, const int pbasis){

  //for(size_t pidx = 0;pidx < pbasis;pidx += nTPB)
  {
  size_t pidx = blockIdx.x*nTPB;
      double sum = 0.0;
      int poffset = pidx + threadIdx.x;

      for(size_t ns = 0;ns < numst;ns++)
      {
          size_t offset = ns*pbasis + poffset;
          if(poffset < pbasis)
              sum += std::real(idata1[offset] * std::conj(idata2[offset]));   
      }
      if(poffset < pbasis)
          odata[poffset] = sum;
  }
}

void GpuProductBr(double *in1, double *in2, double *out, int numst, int pbasis)
{
  int nblocks = pbasis / nTPB + 1;
  breduce<<<nblocks, nTPB>>>(in1, in2, out, numst, pbasis);
}
void GpuProductBr(std::complex<double> *in1, std::complex<double> *in2, double *out, int numst, int pbasis)
{
  int nblocks = pbasis / nTPB + 1;
  breduce<<<nblocks, nTPB>>>(in1, in2, out, numst, pbasis);
}

#endif
