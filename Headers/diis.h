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

#ifndef RMG_diis_H
#define RMG_diis_H 1

#include <vector>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <cmath>
#include <complex>
#include <typeinfo>
#include "GlobalSums.h"
#include "blas.h"
 
template <typename T> class diis {
public:
    diis(int max_Min, int N_in);
    int max_M;
    int N;
    double eps = 1.0e-11;           // epsilon added to B_ii for stability
    void addfunc(T *f);
    void addres(T *r);
    void addres(float *r);
    void addres(std::complex<float> *r);

    std::vector<T> compute_estimate();

    size_t history_size() const { return funcs.size(); }
    void clear() { funcs.clear(); res.clear(); }

private:
    std::vector<std::vector<T>> funcs;
    std::vector<std::vector<T>> res;
    
    T dot(std::vector<T>& a, std::vector<T>& b);

};

#endif
