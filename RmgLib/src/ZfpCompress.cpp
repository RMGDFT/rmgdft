/*
 *
 * Copyright (c) 2018, Emil Briggs
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


#include <typeinfo>
#include "ZfpCompress.h"
template size_t ZfpCompress::compress_buffer<double>(double *, double *, int, int, int, int, size_t);
template size_t ZfpCompress::decompress_buffer<double>(double *, double *, int, int, int, int, size_t);
template size_t ZfpCompress::compress_buffer<float>(float *, float *, int, int, int, int, size_t);
template size_t ZfpCompress::decompress_buffer<float>(float *, float *, int, int, int, int, size_t);

template size_t ZfpCompress::compress_buffer<double>(double *, double *, int, int, int, int, int, int, int, size_t);
template size_t ZfpCompress::decompress_buffer<double>(double *, double *, int, int, int, int, int, int, int, size_t);
template size_t ZfpCompress::compress_buffer<float>(float *, float *, int, int, int, int, int, int, int, size_t);
template size_t ZfpCompress::decompress_buffer<float>(float *, float *, int, int, int, int, int, int, int, size_t);

template size_t ZfpCompress::compress_buffer<double>(double *, double *, int, int, int, double, size_t);
template size_t ZfpCompress::decompress_buffer<double>(double *, double *, int, int, int, double, size_t);
template size_t ZfpCompress::compress_buffer<float>(float *, float *, int, int, int, double, size_t);
template size_t ZfpCompress::decompress_buffer<float>(float *, float *, int, int, int, double, size_t);

template size_t ZfpCompress::compress_buffer<double>(double *, double *, int, int, int, int, int, int, double, size_t);
template size_t ZfpCompress::decompress_buffer<double>(double *, double *, int, int, int, int, int, int, double, size_t);
template size_t ZfpCompress::compress_buffer<float>(float *, float *, int, int, int, int, int, int, double, size_t);
template size_t ZfpCompress::decompress_buffer<float>(float *, float *, int, int, int, int, int, int, double, size_t);

ZfpCompress::ZfpCompress(void)
{
    this->zfp = zfp_stream_open(NULL);
    this->field = zfp_field_alloc();
}

ZfpCompress::~ZfpCompress(void)
{
    zfp_field_free(this->field);
    zfp_stream_close(this->zfp);

}

// Precision variant
template <typename RmgType>
size_t ZfpCompress::compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, zdim, ydim, xdim);
    zfp_field_set_stride_3d(this->field, 1, zdim, ydim*zdim);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_precision(zfp, precision);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Compressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t compressed_size = zfp_compress(this->zfp, this->field);
    if(!compressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error compressing buffer.\n");

    stream_close(stream);
    return compressed_size;
}

// Accuracy variant
template <typename RmgType>
size_t ZfpCompress::compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, double tolerance, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, zdim, ydim, xdim);
    zfp_field_set_stride_3d(this->field, 1, zdim, ydim*zdim);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_accuracy(zfp, tolerance);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Compressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t compressed_size = zfp_compress(this->zfp, this->field);
    if(!compressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error compressing buffer.\n");

    stream_close(stream);
    return compressed_size;
}

// Precision variant
template <typename RmgType>
size_t ZfpCompress::decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, zdim, ydim, xdim);
    zfp_field_set_stride_3d(this->field, 1, zdim, ydim*zdim);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_precision(zfp, precision);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Decompressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t decompressed_size = zfp_decompress(this->zfp, this->field);
    if(!decompressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error decompressing buffer.\n");

    stream_close(stream);

    return decompressed_size;
}

// Accuracy variant
template <typename RmgType>
size_t ZfpCompress::decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, double tolerance, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, zdim, ydim, xdim);
    zfp_field_set_stride_3d(this->field, 1, zdim, ydim*zdim);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_accuracy(zfp, tolerance);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Decompressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t decompressed_size = zfp_decompress(this->zfp, this->field);
    if(!decompressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error decompressing buffer.\n");

    stream_close(stream);

    return decompressed_size;

}



// Strided precision variant
template <typename RmgType>
size_t ZfpCompress::compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, int precision, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, xdim, ydim, zdim);
    zfp_field_set_stride_3d(this->field, sx, sz, sy*sz);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_precision(zfp, precision);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Compressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t compressed_size = zfp_compress(this->zfp, this->field);
    if(!compressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error compressing buffer.\n");

    stream_close(stream);
    return compressed_size;
}

// Strided accuracy variant
template <typename RmgType>
size_t ZfpCompress::compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, double tolerance, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, xdim, ydim, zdim);
    zfp_field_set_stride_3d(this->field, sx, sz, sy*sz);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_accuracy(zfp, tolerance);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Compressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t compressed_size = zfp_compress(this->zfp, this->field);
    if(!compressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error compressing buffer.\n");

    stream_close(stream);
    return compressed_size;
}

// Strided precision variant
template <typename RmgType>
size_t ZfpCompress::decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, int precision, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, xdim, ydim, zdim);
    zfp_field_set_stride_3d(this->field, sx, sz, sy*sz);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_precision(zfp, precision);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Decompressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t decompressed_size = zfp_decompress(this->zfp, this->field);
    if(!decompressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error decompressing buffer.\n");

    stream_close(stream);

    return decompressed_size;

}

// Strided accuracy variant
template <typename RmgType>
size_t ZfpCompress::decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, double tolerance, size_t outbufsize)
{
    if(typeid(RmgType) == typeid(double)) zfp_field_set_type(this->field, zfp_type_double);
    if(typeid(RmgType) == typeid(float)) zfp_field_set_type(this->field, zfp_type_float);

    zfp_field_set_size_3d(this->field, xdim, ydim, zdim);
    zfp_field_set_stride_3d(this->field, sx, sz, sy*sz);
    zfp_field_set_pointer(this->field, in);
    zfp_stream_set_accuracy(zfp, tolerance);

    size_t bufsize = zfp_stream_maximum_size(this->zfp, this->field);
    if(bufsize > outbufsize)
        rmg_error_handler (__FILE__, __LINE__, "Decompressed buffer is too small.\n");

    bitstream *stream = stream_open(out, bufsize);
    zfp_stream_set_bit_stream(this->zfp, stream);
    zfp_stream_rewind(this->zfp);
    size_t decompressed_size = zfp_decompress(this->zfp, this->field);
    if(!decompressed_size)
        rmg_error_handler (__FILE__, __LINE__, "Error decompressing buffer.\n");

    stream_close(stream);

    return decompressed_size;

}
