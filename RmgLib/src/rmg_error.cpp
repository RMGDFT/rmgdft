/*
 *
 * Copyright (c) 2014, Emil Briggs
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


#include "rmg_error.h"
#include <iostream>
#include <signal.h>

#if HIP_ENABLED
#include "tiled_mm.hpp"
#endif  

#include <boost/stacktrace.hpp>

static int do_print=true;

void rmg::error_set_print(int doprint)
{
    do_print = doprint;
}

void rmg::error(char const *message, const std::source_location loc)
{
   if(do_print)
   {
        std::cout << "\nError: " << message << "\n";
        std::cout << "Function:  " << loc.function_name() << "\n"
                  << "File:      " << loc.file_name()     << "\n"
                  << "Line:      " << loc.line()          << "\n\n";
        std::cout << "Stacktrace:" << "\n";
        std::cout << boost::stacktrace::stacktrace();
    }
    fflush (NULL);
    sleep(1);
    raise(SIGTERM);
}

#if CUDA_ENABLED
void rmg::error(cublasStatus_t custat, const std::source_location loc)
{

    if(custat==CUBLAS_STATUS_SUCCESS)
        return;

    std::string msg;
    if(custat==CUBLAS_STATUS_NOT_INITIALIZED)
    {
        msg = "CUBLAS_STATUS_NOT_INITIALIZED'";
    }
    else if(custat==CUBLAS_STATUS_ALLOC_FAILED)
    {
        msg = "CUBLAS_STATUS_ALLOC_FAILED";
    }
    else if(custat==CUBLAS_STATUS_INVALID_VALUE)
    {
        msg = "CUBLAS_STATUS_INVALID_VALUE";
    }
    else if(custat==CUBLAS_STATUS_ARCH_MISMATCH)
    {
        msg = "CUBLAS_STATUS_ARCH_MISMATCH";
    }
    else if(custat==CUBLAS_STATUS_MAPPING_ERROR)
    {
        msg = "CUBLAS_STATUS_MAPPING_ERROR";
    }
    else if(custat==CUBLAS_STATUS_EXECUTION_FAILED)
    {
        msg = "CUBLAS_STATUS_EXECUTION_FAILED";
    }
    else if(custat==CUBLAS_STATUS_INTERNAL_ERROR)
    {
        msg = "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    else
    {
        msg = "UNKNOWN CUBLAS ERROR";
    }

    if(do_print)
    {
        std::cout << "\nError: " << msg << "\n";
        std::cout << "Function:  " << loc.function_name() << "\n"
                  << "File:      " << loc.file_name()     << "\n"
                  << "Line:      " << loc.line()          << "\n\n";
        std::cout << "Stacktrace:" << "\n";
        std::cout << boost::stacktrace::stacktrace();
    }
    fflush (NULL);
    sleep(1);
    raise(SIGTERM);

}

void rmg::error(cudaError_t custat, const std::source_location loc)
{

    if(custat == cudaSuccess) return;

    const char *msg = cudaGetErrorString(custat);
    if(do_print)
    {
        std::cout << "\nError: " << msg << "\n";
        std::cout << "Function:  " << loc.function_name() << "\n"
                  << "File:      " << loc.file_name()     << "\n"
                  << "Line:      " << loc.line()          << "\n\n";
        std::cout << "Stacktrace:" << "\n";
        std::cout << boost::stacktrace::stacktrace();
    }
    fflush (NULL);
    sleep(1);
    raise(SIGTERM);
}

#endif

#if HIP_ENABLED
void rmg::error(hipblasStatus_t hipstat, const std::source_location loc)
{

    if(hiptat==HIPBLAS_STATUS_SUCCESS)
        return;

    std::string msg;
    if(hipstat==HIPBLAS_STATUS_NOT_INITIALIZED)
    {
        msg = "'HIPBLAS_STATUS_NOT_INITIALIZED'";
    }
    else if(hipstat==HIPBLAS_STATUS_ALLOC_FAILED)
    {
        msg = "HIPBLAS_STATUS_ALLOC_FAILED";
    }
    else if(hipstat==HIPBLAS_STATUS_INVALID_VALUE)
    {
        msg = "HIPBLAS_STATUS_INVALID_VALUE";
    }
    else if(hipstat==HIPBLAS_STATUS_ARCH_MISMATCH)
    {
        msg = "HIPBLAS_STATUS_ARCH_MISMATCH";
    }
    else if(hipstat==HIPBLAS_STATUS_MAPPING_ERROR)
    {
        msg = "HIPBLAS_STATUS_MAPPING_ERROR";
    }
    else if(hipstat==HIPBLAS_STATUS_EXECUTION_FAILED)
    {
        msg = "HIPBLAS_STATUS_EXECUTION_FAILED";
    }
    else if(hipstat==HIPBLAS_STATUS_INTERNAL_ERROR)
    {
        msg = "HIPBLAS_STATUS_INTERNAL_ERROR";
    }
    else
    {
        msg = "UNKNOWN HIPBLAS ERROR";
    }

    if(do_print)
    {
        std::cout << "\nError: " << msg << "\n";
        std::cout << "Function:  " << loc.function_name() << "\n"
                  << "File:      " << loc.file_name()     << "\n"
                  << "Line:      " << loc.line()          << "\n\n";
        std::cout << "Stacktrace:" << "\n";
        std::cout << boost::stacktrace::stacktrace();
    }
    fflush (NULL);
    sleep(1);
    raise(SIGTERM);

}

void rmg::error(hipError_t hipstat, const std::source_location loc)
{

    if(hipstat == hipSuccess) return;

    const char *msg = hipGetErrorString(hipstat);
    if(do_print)
    {   
        std::cout << "\nError: " << msg << "\n";
        std::cout << "Function:  " << loc.function_name() << "\n"
                  << "File:      " << loc.file_name()     << "\n"
                  << "Line:      " << loc.line()          << "\n\n";
        std::cout << "Stacktrace:" << "\n";
        std::cout << boost::stacktrace::stacktrace();
    }
    fflush (NULL);
    sleep(1);
    raise(SIGTERM);
}   
        


#endif
