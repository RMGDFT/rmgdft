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

#ifndef RMG_ZfpCompress_H
#define RMG_ZfpCompress_H 1

#include "rmg_error.h"
#include "zfp.h"

#ifdef __cplusplus


class ZfpCompress {

protected:
  
private:
    zfp_stream *zfp;   /* compressed stream */
    zfp_field *field;  /* array meta data */

public:

    ZfpCompress(void);
    ~ZfpCompress(void);

    template <typename RmgType> size_t compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize);
    template <typename RmgType> size_t decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize);
    template <typename RmgType> size_t compress_buffer_fixed_rate(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize);
    template <typename RmgType> size_t decompress_buffer_fixed_rate(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int precision, size_t outbufsize);
    template <typename RmgType> size_t compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, double tolerance, size_t outbufsize);
    template <typename RmgType> size_t decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, double tolerance, size_t outbufsize);

    template <typename RmgType> size_t compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, int precision, size_t outbufsize);
    template <typename RmgType> size_t decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, int precision, size_t outbufsize);
    template <typename RmgType> size_t compress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, double tolerance, size_t outbufsize);
    template <typename RmgType> size_t decompress_buffer(RmgType *in, RmgType *out, int xdim, int ydim, int zdim, int sx, int sy, int sz, double tolerance, size_t outbufsize);
};

#endif
#endif
