/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
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



#if CUDA_ENABLED
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include "Gpufuncs.h"

// This is used in the Folded Spectrum routine. If an element of dx
// which is a pointer to a block of memory allocated by GpuMallocManged
// is less than zero it reverses the sign of the corresponding element
// in dy.
template <typename Iterator>
class strided_range
{
    public:

    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

    struct stride_functor : public thrust::unary_function<difference_type,difference_type>
    {
        difference_type stride;

        stride_functor(difference_type stride)
            : stride(stride) {}

        __host__ __device__
        difference_type operator()(const difference_type& i) const
        {
            return stride * i;
        }
    };

    typedef typename thrust::counting_iterator<difference_type>                   CountingIterator;
    typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
    typedef typename thrust::permutation_iterator<Iterator,TransformIterator>     PermutationIterator;

    // type of the strided_range iterator
    typedef PermutationIterator iterator;

    // construct strided_range for the range [first,last)
    strided_range(Iterator first, Iterator last, difference_type stride)
        : first(first), last(last), stride(stride) {}

    iterator begin(void) const
    {
        return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
    }

    iterator end(void) const
    {
        return begin() + ((last - first) + (stride - 1)) / stride;
    }

    protected:
    Iterator first;
    Iterator last;
    difference_type stride;
};
struct is_negative
{
  __host__ __device__
  bool operator()(double x)
  {
    return (x < 0.0);
  }
};

__global__ void NegateDiag(double *dx, int incx, double *dy, int incy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) 
    {
        if(dx[i*incx + i] < 0.0) dy[i*incy] = -dy[i*incy];
    }
}


void GpuNegate(double *dx, int incx, double *dy, int incy, int n)
{
#if 0
    DeviceSynchronize();
    thrust::negate<double> neg_op;
    thrust::device_ptr<double> dxptr = thrust::device_pointer_cast(dx);
    thrust::device_ptr<double> dyptr = thrust::device_pointer_cast(dy);
    typedef thrust::device_vector<double>::iterator Iterator;
    strided_range<Iterator> pos(dxptr, dxptr + n, incx);
    //thrust::transform_if(dxptr, dxptr + n, dyptr, neg_op, is_negative());
    thrust::transform_if(pos.begin(), pos.end(), dyptr, neg_op, is_negative());
    DeviceSynchronize();
#else
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / n;
    NegateDiag<<<numBlocks, blockSize>>>(dx, incx, dy, incy, n);
#endif

}

#endif
