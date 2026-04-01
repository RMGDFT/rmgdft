/*
 *
 * Copyright (c) 2026, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/


#include "rmg_hvector.h"

namespace rmg
{
template<typename T>
hvector<T>::hvector(size_t pbasis_in)
{
    pbasis_ = pbasis_in;
#if CUDA_ENABLED || HIP_ENABLED
    hpool.malloc(&data_, pbasis_in);
#elif SYCL_ENABLED
    data_ = sycl::malloc_host<T>(pbasis_in, hpool.sycl_Q);
#else
    data_ = new T[pbasis_in];
#endif
}

template<typename T>
hvector<T>::~hvector(void)
{
#if CUDA_ENABLED || HIP_ENABLED
    hpool.free(data_);
#elif SYCL_ENABLED
    sycl::free(data_, hpool.sycl_Q);
#else
    delete [] data_;
#endif
}


// Instantiate all versions
template hvector<float>::hvector(size_t);
template hvector<double>::hvector(size_t);
template hvector<std::complex<float>>::hvector(size_t);
template hvector<std::complex<double>>::hvector(size_t);
template hvector<float>::~hvector(void);
template hvector<double>::~hvector(void);
template hvector<std::complex<float>>::~hvector(void);
template hvector<std::complex<double>>::~hvector(void);

template void hvector<float>::multiply(const float&);
template void hvector<double>::multiply(const double&);
template void hvector<std::complex<float>>::multiply(const std::complex<float>&);
template void hvector<std::complex<double>>::multiply(const std::complex<double>&);
}
