/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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

#include "rmg_tddft.h"
#include "../Headers/prototypes_tddft.h"
#include "rmg_dev_allocate.h"
#include "Kpoint.h"
#include "blas_driver.h"


template void rmg::rotate_sint<double, double>(Kpoint<double> *, double *, double *, int, int);
template void rmg::rotate_sint<double, std::complex<double>>(Kpoint<double> *, double *, std::complex<double> *, int, int);
template void rmg::rotate_sint<std::complex<double>, std::complex<double>>(Kpoint<std::complex<double>> *, std::complex<double> *, std::complex<double> *, int, int);
template <typename OrbitalType, typename MatrixType>

void rmg::rotate_sint(
            Kpoint<OrbitalType> *Kptr, OrbitalType *sint, MatrixType *rho_matrix, int offset, int nstates)
{
}
