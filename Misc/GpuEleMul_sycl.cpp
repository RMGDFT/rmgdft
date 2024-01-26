/*
 *
 * Copyright 2023 The RMG Project Developers. See the COPYRIGHT file 
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



#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include <sycl/queue.hpp>
    #include "Gpufuncs.h"


void GpuEleMul(double *dx, std::complex<double> *dy, int n, cl::sycl::queue &Q)
{
    sycl::range<1> num_items{(size_t)n};
    auto e= Q.parallel_for(num_items, [=](auto i) { dy[i] = dx[i] * dy[i]; });
    e.wait();
}

void GpuEleMul(double *dx, std::complex<float> *dy, int n, cl::sycl::queue &Q)
{
    sycl::range<1> num_items{(size_t)n};
    auto e= Q.parallel_for(num_items, [=](auto i) { dy[i] = ((float)dx[i]) * dy[i]; });
    e.wait();
}

#endif
