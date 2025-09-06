/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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

#include <vector>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <cmath>
#include <complex>
#include <typeinfo>
#include "transition.h"
#include "const.h"
#include "pe_control.h"
#include "rmg_control.h"
#include "Solvers.h"



// We use this routine instead of lapack in cases where the matrices are very small
// and performance is not an issue but the function is called from a threaded region
// as differences between lapack and blas implementations may lead to problems.
template <typename T> void gauss_solve(std::vector<T>& A, std::vector<T>& b, int N)
{
    const double eps = 1.0e-15;

    for (int k = 0; k < N; ++k) {
        int ipiv = k;
        double maxabs = std::abs(A[k * N + k]);
        for (int i = k + 1; i < N; ++i)
        {
            double v = std::abs(A[i * N + k]);
            if (v > maxabs) { maxabs = v; ipiv = i; }
        }
        if (maxabs < eps)
        {
            rmg_error_handler (__FILE__,__LINE__,"Singular matrix encountered.");
        }
        if (ipiv != k)
        {
            for (int j = k; j < N; ++j) std::swap(A[k * N + j], A[ipiv * N + j]);
            std::swap(b[k], b[ipiv]);
        }

        T Akk = A[k * N + k];
        for (int i = k + 1; i < N; ++i)
        {
            T f = A[i * N + k] / Akk;
            if (std::abs(f) < eps) continue;
            A[i * N + k] = 0.0;
            for (int j = k + 1; j < N; ++j) A[i * N + j] -= f * A[k * N + j];
            b[i] -= f * b[k];
        }
    }

    for (int i = N - 1; i >= 0; --i)
    {
        T s = b[i];
        for (int j = i + 1; j < N; ++j)
            s -= A[i * N + j] * b[j];
        T Aii = A[i * N + i];
        if (std::abs(Aii) < eps)
        {
            rmg_error_handler (__FILE__,__LINE__,"Singular matrix encountered.");
        }
        b[i] = s / Aii;
    }
}

template void gauss_solve<double>(std::vector<double>& A, std::vector<double>& b, int N);
template void gauss_solve<std::complex<double>>(std::vector<std::complex<double>>& A, std::vector<std::complex<double>>& b, int N);

