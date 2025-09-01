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
#include "GlobalSums.h"
#include "pe_control.h"
#include "rmg_control.h"
#include "blas.h"
#include "diis.h"

template <class T> diis<T>::diis(int max_Min, int N_in)
{
    max_M = max_Min;      // Maximum number of past iterates to use
    N = N_in;             // array length for functions and residuals
}


// Adds an iteration of the function to be optimized.
template <class T> void diis<T>::addfunc(T *f)
{
    std::vector<T> ftmp(N);
    for(int i=0;i < N;i++) ftmp[i] = f[i];
    funcs.push_back(ftmp);

    // Remove oldest entry if needed
    if (funcs.size() > max_M) {
        funcs.erase(funcs.begin());
    }
}

template <class T> void diis<T>::addfunc(float *f)
{
    if constexpr (std::is_same_v<T, double>)
    {
        std::vector<T> ftmp(N);
        for(int i=0;i < N;i++) ftmp[i] = f[i];
        funcs.push_back(ftmp);

        // Remove oldest entry if needed
        if (funcs.size() > max_M) {
            funcs.erase(funcs.begin());
        }
    }
}

template <class T> void diis<T>::addfunc(std::complex<float> *f)
{
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        std::vector<T> ftmp(N);
        for(int i=0;i < N;i++) ftmp[i] = f[i];
        funcs.push_back(ftmp);

        // Remove oldest entry if needed
        if (funcs.size() > max_M) {
            funcs.erase(funcs.begin());
        }
    }
}

// Adds a double or std::complex<double> residual
// Adds a double or std::complex<double> residual
template <class T> void diis<T>::addres(T *r)
{
    if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
    {
        std::vector<T> rtmp(N);
        for(int i=0;i < N;i++) rtmp[i] = r[i];
        res.push_back(std::move(rtmp));

        // Remove oldest entry if needed
        if (res.size() > max_M) {
            res.erase(res.begin());
        }
    }
}

// Adds a float residual
template <class T> void diis<T>::addres(float *r)
{
    if constexpr (std::is_same_v<T, double>)
    {
        std::vector<T> rtmp(N);
        for(int i=0;i < N;i++) rtmp[i] = r[i];
        res.push_back(std::move(rtmp));

        // Remove oldest entry if needed
        if (res.size() > max_M) {
            res.erase(res.begin());
        }
    } 
}

// Adds a std::complex<float> residual
template <class T> void diis<T>::addres(std::complex<float> *r)
{
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        std::vector<T> rtmp(N);
        for(int i=0;i < N;i++) rtmp[i] = r[i];
        res.push_back(std::move(rtmp));

        // Remove oldest entry if needed
        if (res.size() > max_M) {
            res.erase(res.begin());
        }
    } 
}

// Compute the DIIS estimated state. Requires at least 2 history entries.
// Returns the mixed state. Throws if history is insufficient.
template <class T> std::vector<T> diis<T>::compute_estimate()
{
    int m = res.size();
    if (m < 2)
        rmg_error_handler(__FILE__, __LINE__,"diis requires at least 2 entries.");

    // Generate square matrix B(M,M)
    // B_ij = <r_i, r_j> for ij < M-1
    // and setup coeffs so that c_i sum to 1.0
    int M = m + 1;
    std::vector<T> A(M * M, 0.0);
    std::vector<T> b(M, 0.0);
    b[m] = 1.0;
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    T vel = get_vel();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= i; j++) {
            //T val = dot(res[i], res[j]);
           T val = 0.0;
           if constexpr (std::is_same_v<T, double>)
               for (int k = 0; k < N; ++k) val = val + res[i][k] * res[j][k];
           else
               for (int k = 0; k < N; ++k) val = val + std::conj(res[i][k]) * res[j][k];

            GlobalSums(&val, factor, pct.coalesced_grid_comm);
            A[i * M + j] = val;
            A[j * M + i] = val;
        }
    }
//        GlobalSums((double *)A.data(), factor*m*m, pct.grid_comm);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            A[j * M + i] *= vel;
        }
    }

    // Find maxdiag as a proxy for a singular matrix
    // and add eps to diagonal entries except the last
    double maxdiag = 0.0;
    for (int i = 0; i < m; ++i) {
        maxdiag = std::max(maxdiag, std::abs(A[i*M+i]));
// May be better to just bail if diagonal entries are smaller than eps
        A[i*M + i] += eps;
        A[i * M + m] = 1.0;
        A[m * M + i] = 1.0;
    }
    A[m * M + m] = 0.0;

    // bail if matrix becomes singular?
    // use zero length return vector if so
//    eps = std::min(1.0e-12, ct.scf_accuracy/100000.0);
    if(cleared || (maxdiag < eps)) {
        std::vector<T> rm;
        rm.clear();
        cleared = true;
        //if(pct.gridpe==0)printf("CLEARED 0\n");
        return rm;
    }

    // Solve A*x = b with lapack
    std::vector<T> x = b;
    std::vector<int> ipvt(N);
    int ione = 1, info;
    if(typeid(T) == typeid(double))
    {
#if 0
        dgesv(&M, &ione, (double *)A.data(), &M, ipvt.data(), (double *)b.data(), &M, &info);
#else
    // dgelss better for nearly singular matrices?
    int ione = 1, rank, info;
    int lwork=100;
    double work[100];
    std::vector<double> S(M);    // singular values
    double rcond = -1.0;
    dgelss(&M, &M, &ione, (double *)A.data(), &M,
                     (double *)b.data(), &M,
                     (double *)S.data(), &rcond,
                     &rank,
                     (double *)work, &lwork,
                     &info);
#endif
    }
    else
    {
        zgesv(&M, &ione, (double *)A.data(), &M, ipvt.data(),
             (double *)b.data(), &M, &info);
    }

    // Form the new estimate using the coefficients c_i
    std::vector<T> mixed = funcs.back();
    std::fill(mixed.begin(), mixed.end(), 0.0);
    double maxci = 0.0;
    for (int i = 0; i < m; i++)
    {
        T c = b[i];
        maxci = std::max(maxci, std::abs(c));
        const auto& t = funcs[i];
        for (int k = 0; k < N; k++) mixed[k] += c * t[k];
    }
    if(cleared || maxci > 2.0*m)
    {
        std::vector<T> rm;
        rm.clear();
        cleared = true;
        //if(pct.gridpe==0)printf("CLEARED 1\n");
        return rm;
    }

    // Make sure the new wavefunction is normalized
    double n = std::sqrt(get_vel()*std::real(dot(mixed, mixed)));
    double invm = 1.0/(n + 1.0e-30);
    for (auto& v : mixed) v *= invm;
    return mixed;
}

template <class T> T diis<T>::dot(std::vector<T>& a, std::vector<T>& b)
{
    int factor = 1;
    if(typeid(T) == typeid(std::complex<double>)) factor = 2;
    T sum = 0.0;
    if constexpr (std::is_same_v<T, std::complex<double>>)
        for (int i = 0; i < N; ++i) sum += std::conj(a[i]) * b[i];
    else
        for (int i = 0; i < N; ++i) sum += a[i] * b[i];
    GlobalSums((double *)&sum, factor, pct.coalesced_grid_comm);
    return sum;
}


template diis<double>::diis(int max_Nin, int N_in);
template diis<std::complex<double>>::diis(int max_Nin, int N_in);
template void diis<double>::addfunc(double *f);
template void diis<double>::addfunc(float *f);
template void diis<std::complex<double>>::addfunc(std::complex<double> *f);
template void diis<std::complex<double>>::addfunc(std::complex<float> *f);
template void diis<double>::addres(double *f);
template void diis<double>::addres(float *f);
template void diis<std::complex<double>>::addres(std::complex<double> *f);
template void diis<std::complex<double>>::addres(std::complex<float> *f);
template std::vector<double> diis<double>::compute_estimate(void);
template std::vector<std::complex<double>> diis<std::complex<double>>::compute_estimate(void);

