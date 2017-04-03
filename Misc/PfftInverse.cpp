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
#include "Pw.h"
#include "transition.h"
#include "RmgParallelFft.h"


// This function performs a parallel inverse FFT of a std::complex<double> array using the pfft code from lammps
void PfftInverse (std::complex<double> * in, std::complex<double> * out, Pw &pwaves)
{

  std::complex<double> *buf = new std::complex<double>[pwaves.pbasis];
  for(int i = 0;i < pwaves.pbasis;i++) buf[i] = in[i];
  fft_3d((FFT_DATA *)buf, (FFT_DATA *)out, 1, pwaves.fft_backward_plan);

  delete [] buf;
}


