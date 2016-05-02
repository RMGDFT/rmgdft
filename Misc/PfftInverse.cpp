/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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
#include <math.h>
#include <float.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "Pw.h"
#include "transition.h"
#include "RmgParallelFft.h"


// This function performs a parallel inverse FFT of a std::complex<double> array using the pfft library.
// If the actual processor grid is not directly compatible with the pfft requirements
// it handles the required data movement in order to perform the FFT. No scaling is done.
void PfftInverse (std::complex<double> * in, std::complex<double> * out, Pw &pwaves)
{

  int size = std::max(pwaves.pbasis, pwaves.remap_local_size);
  std::complex<double> *buf = new std::complex<double>[size];
  for(int i = 0;i < pwaves.pbasis;i++) buf[i] = in[i];

  if(pwaves.fwd_remap) remap_3d((double *)buf, (double *)buf, NULL, pwaves.fwd_remap);
  pfft_execute_dft(*pwaves.backward_plan, (double (*)[2])buf, (double (*)[2])buf);
  if(pwaves.inv_remap) remap_3d((double *)buf, (double *)buf, NULL, pwaves.inv_remap);

  for(int i = 0;i < pwaves.pbasis;i++) out[i] = buf[i];

  delete [] buf;
}




