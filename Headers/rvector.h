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


/* This class implements rvectors for RMG. An rvector is a memory block that is allocated using an upstream
   memory resource when the rvector is created. It is automatically deallocated when the object goes out
   of scope.


   Constructor examples:
     rvector<T> b(size);
     Creates an revector object b which uses a memory allocation of size*sizeof(T) bytes.
     the object is deleted or goes out of scope.

*/



#ifndef RMG_rvector_H
#define RMG_rvector_H 1

#include <complex>
#include <span>
#include "rmg_allocate.h"


template <typename T> class rvector {

  friend rvector operator+(rvector a, const rvector& b) {
    a += b;
    return std::move(a);
  }
  friend rvector& operator+=(rvector& a, const rvector& b) {
    a.increment(b);
    return a;
  }
  friend rvector operator-(rvector a, const rvector& b) {
    a -= b;
    return std::move(a);
  }
  friend rvector& operator-=(rvector& a, const rvector& b) {
    a.decrement(b);
    return a;
  }
  friend rvector operator*(const T &b, rvector a) {
    a.multiply(b);
    return a;
  }
  friend rvector operator*(rvector a, const T &b) {
    a *= b;
    return a;
  }
  friend rvector& operator*=(rvector& a, const T &b) {
    a.multiply(b);
    return a;
  }

public:
    rvector(size_t pbasis_in);
    rvector(rvector& t) // copy constructor
    {
        this->pbasis_   = t.pbasis;
        this->allocate();
        std::copy(t.data_, t.data_ + this->pbasis_, this->data_);
    }
    ~rvector(void);

    void set(T x)
    {
        for(int i=0;i < this->pbasis_;i++) this->data_[i] = x;
    }

    const int& pbasis = pbasis_;

    const int size() const { return pbasis; }
    T* data() { return data_; }

    T& operator [](int idx) {
        return data_[idx];
    }

    T operator [](int idx) const {
        return data_[idx];
    }

private:

protected:
   int pbasis_;
   T *data_;

   void allocate(void)
   {
#if CUDA_ENABLED || HIP_ENABLED
#else
       data_ = new T[pbasis_]();
#endif

   }

void increment(const rvector<T>& c) {
  if(this->pbasis == c.pbasis) {
    for (int i = 0; i < this->pbasis; i++) {
      data_[i] += c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

void decrement(const rvector<T>& c) {
  if(this->pbasis == c.pbasis) {
    for (int i = 0; i < this->pbasis; i++) {
      data_[i] -= c.data_[i];
    }
  } else {
    throw "rvector objects are not the same size!";
  }
}

void multiply(const T& b) {
    for (int i = 0; i < this->pbasis; i++) {
        data_[i] *= b;
    }
}

   rvector& operator=(rvector const &rhs)
   {
       if (this != &rhs)
       {
           std::copy(rhs.data_, rhs.data_ + this->pbasis_, this->data_);
       }
       return *this;
    };

};


#endif
