/*
 *
 * Copyright (c) 2024, Emil Briggs
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

#include "GridObject.h"

template<typename T>
GridObject<T>::GridObject(int density)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    pbasis = dimx_ * dimy_ * dimz_;
    data_ = new T[pbasis]();
}

template<typename T>
GridObject<T>::GridObject(int density, T *data_ptr)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    pbasis = dimx_ * dimy_ * dimz_;
    data_ = data_ptr;
    owns_allocation = false;
}

template<typename T>
GridObject<T>::~GridObject(void)
{
    if(owns_allocation) delete [] data_;
}

template<typename T>
void GridObject<T>::increment(const GridObject<T>& c) {
  if(this->pbasis == c.pbasis) {
    for (int i = 0; i < this->pbasis; i++) {
      data_[i] += c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::decrement(const GridObject<T>& c) {
  if(this->pbasis == c.pbasis) {
    for (int i = 0; i < this->pbasis; i++) {
      data_[i] -= c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::multiply(const T& b) {
    for (int i = 0; i < this->pbasis; i++) {
        data_[i] *= b;
    }
}

// Instantiate all versions
template GridObject<float>::GridObject(int);
template GridObject<double>::GridObject(int);
template GridObject<std::complex<float>>::GridObject(int);
template GridObject<std::complex<double>>::GridObject(int);
template GridObject<float>::~GridObject(void);
template GridObject<double>::~GridObject(void);
template GridObject<std::complex<float>>::~GridObject(void);
template GridObject<std::complex<double>>::~GridObject(void);

template FineGridObject<float>::FineGridObject(void);
template FineGridObject<double>::FineGridObject(void);
template FineGridObject<std::complex<float>>::FineGridObject(void);
template FineGridObject<std::complex<double>>::FineGridObject(void);
template FineGridObject<float>::~FineGridObject(void);
template FineGridObject<double>::~FineGridObject(void);
template FineGridObject<std::complex<float>>::~FineGridObject(void);
template FineGridObject<std::complex<double>>::~FineGridObject(void);

