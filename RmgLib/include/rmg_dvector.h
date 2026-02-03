#pragma once
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


/*





*/


/* This class implements dvectors for RMG. A dvector is a memory block that resides on a GPU device
   and is allocated from a pool. It is automatically deallocated when the object goes out
   of scope.

   Constructor examples:
     dvector<T> b(size);
     Creates an dvector object b which uses a memory allocation of size*sizeof(T) bytes.

*/



#ifndef RMG_dvector_H
#define RMG_dvector_H 1

#include <complex>
#include <span>
#include "rmg_dev_allocate.h"

namespace rmg
{
    template <typename T> class dvector {

    public:
        dvector(size_t pbasis_in);
#if 0
        dvector(dvector& t) // copy constructor
        {
            this->pbasis_   = t.pbasis;
            this->allocate();
            std::copy(t.data_, t.data_ + this->pbasis_, this->data_);
        }
#endif
        ~dvector(void);

        void set(T x)
        {
         //   for(int i=0;i < this->pbasis_;i++) this->data_[i] = x;
        }

        const int& pbasis = pbasis_;

        const int size() const { return pbasis; }
        T* data() { return data_; }

    private:

    protected:
       int pbasis_;
       T *data_;
    };
}

#endif
