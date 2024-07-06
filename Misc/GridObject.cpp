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

#include "transition.h"
#include "GridObject.h"

template<typename T>
GridObject<T>::GridObject(int density)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    incz_ = 1;
    incy_ = dimz_;
    incx_ = dimy_*dimz_;
    pbasis_ = dimx_ * dimy_ * dimz_;
}

template<typename T>
GridObject<T>::GridObject(int density, T *data_ptr)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    incz_ = 1;
    incy_ = dimz_;
    incx_ = dimy_*dimz_;
    offsetx_ = Rmg_G->get_PX_OFFSET(density);
    offsety_ = Rmg_G->get_PY_OFFSET(density);
    offsetz_ = Rmg_G->get_PZ_OFFSET(density);
    pbasis_ = dimx_ * dimy_ * dimz_;
}

template<typename T>
GridObject<T>::~GridObject(void)
{
    if(owns_allocation) delete [] data_;
}

template<typename T>
void GridObject<T>::increment(const GridObject<T>& c) {
  if(this->pbasis == c.pbasis && this->factor == c.factor) {
    for (int i = 0; i < this->factor*this->pbasis; i++) {
      data_[i] += c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::increment(const fgobj<T>& c) {
  if(this->pbasis == c.pbasis && this->factor == c.factor) {
    for (int i = 0; i < this->factor*this->pbasis; i++) {
      data_[i] += c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::decrement(const GridObject<T>& c) {
  if(this->pbasis == c.pbasis && this->factor == c.factor) {
    for (int i = 0; i < factor*this->pbasis; i++) {
      data_[i] -= c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::decrement(const fgobj<T>& c) {
  if(this->pbasis == c.pbasis && this->factor == c.factor) {
    for (int i = 0; i < factor*this->pbasis; i++) {
      data_[i] -= c.data_[i];
    }
  } else {
    throw "Grid objects are not the same size!";
  }
}

template<typename T>
void GridObject<T>::multiply(const T& b) {
    for (int i = 0; i < factor * this->pbasis; i++) {
        data_[i] *= b;
    }
}

template<typename T>
wfobj<T>::wfobj(void) : GridObject<T>(1)
{   
    // For wfobj factor is 1 or 2 so for 1 up=dw
    this->allocate(ct.noncoll_factor);
    this->up = std::span(this->data_, this->pbasis);
    this->dw = std::span(this->data_ + (this->factor-1)*this->pbasis, this->pbasis);
}

template<typename T>
wfobj<T>::wfobj(T *data_ptr) : GridObject<T>(1)
{   
    // For wfobj factor is 1 or 2 so for 1 up=dw
    this->allocate(ct.noncoll_factor, data_ptr);
    this->up = std::span(this->data_, this->pbasis);
    this->dw = std::span(this->data_ + (this->factor-1)*this->pbasis, this->pbasis);
}

template<typename T>
wfobj<T>::~wfobj(void)
{
}

template<typename T>
fgobj<T>::fgobj(void) : GridObject<T>(Rmg_G->default_FG_RATIO)
{   
    this->allocate(1);
}

template<typename T>
fgobj<T>::fgobj(T *data_ptr) : GridObject<T>(Rmg_G->default_FG_RATIO)
{   
    this->allocate(1, data_ptr);
}

template<typename T>
fgobj<T>::~fgobj(void)
{
}

template<typename T>
spinobj<T>::spinobj(void) : GridObject<T>(Rmg_G->default_FG_RATIO)
{
    this->allocate(ct.nspin);
    up = std::span(this->data_, this->pbasis);
    dw = std::span(this->data_ + (this->factor-1)*this->pbasis, this->pbasis);
    // We don't define the 4 component names for the spin polarized case. This
    // will cause an error if someone tries to use them like this which is fine.
    if(this->factor == 4)
    {
        c0 = std::span(this->data_, this->pbasis);
        cx = std::span(this->data_ + this->pbasis, this->pbasis);
        cy = std::span(this->data_ + 2*this->pbasis, this->pbasis);
        cz = std::span(this->data_ + 3*this->pbasis, this->pbasis);
    }
}
template<typename T>
spinobj<T>::spinobj(T *data_ptr) : GridObject<T>(Rmg_G->default_FG_RATIO)
{
    this->allocate(ct.nspin, data_ptr);
    up = std::span(this->data_, this->pbasis);
    dw = std::span(this->data_ + (this->factor-1)*this->pbasis, this->pbasis);
}
template<typename T>
spinobj<T>::~spinobj(void)
{
}

template<typename T>
void spinobj<T>::get_oppo(void)
{
    if(this->factor != 2) return;
    get_rho_oppo(this->up.data(), this->dw.data());
}

// Instantiate all versions
template GridObject<float>::GridObject(int);
template GridObject<double>::GridObject(int);
template GridObject<float>::GridObject(int, float *);
template GridObject<double>::GridObject(int, double *);
template GridObject<std::complex<float>>::GridObject(int);
template GridObject<std::complex<double>>::GridObject(int);
template GridObject<std::complex<float>>::GridObject(int, std::complex<float> *);
template GridObject<std::complex<double>>::GridObject(int, std::complex<double> *);
template GridObject<float>::~GridObject(void);
template GridObject<double>::~GridObject(void);
template GridObject<std::complex<float>>::~GridObject(void);
template GridObject<std::complex<double>>::~GridObject(void);
template void GridObject<float>::increment(const fgobj<float>&);
template void GridObject<double>::increment(const fgobj<double>&);
template void GridObject<std::complex<float>>::increment(const fgobj<std::complex<float>>&);
template void GridObject<std::complex<double>>::increment(const fgobj<std::complex<double>>&);
template void GridObject<float>::decrement(const fgobj<float>&);
template void GridObject<double>::decrement(const fgobj<double>&);
template void GridObject<std::complex<float>>::decrement(const fgobj<std::complex<float>>&);
template void GridObject<std::complex<double>>::decrement(const fgobj<std::complex<double>>&);

template void GridObject<float>::multiply(const float&);
template void GridObject<double>::multiply(const double&);
template void GridObject<std::complex<float>>::multiply(const std::complex<float>&);
template void GridObject<std::complex<double>>::multiply(const std::complex<double>&);

template fgobj<float>::fgobj(void);
template fgobj<double>::fgobj(void);
template fgobj<std::complex<float>>::fgobj(void);
template fgobj<std::complex<double>>::fgobj(void);
template fgobj<float>::fgobj(float *);
template fgobj<double>::fgobj(double *);
template fgobj<std::complex<float>>::fgobj(std::complex<float> *);
template fgobj<std::complex<double>>::fgobj(std::complex<double> *);
template fgobj<float>::~fgobj(void);
template fgobj<double>::~fgobj(void);
template fgobj<std::complex<float>>::~fgobj(void);
template fgobj<std::complex<double>>::~fgobj(void);

template wfobj<float>::wfobj(void);
template wfobj<double>::wfobj(void);
template wfobj<std::complex<float>>::wfobj(void);
template wfobj<std::complex<double>>::wfobj(void);
template wfobj<float>::wfobj(float *);
template wfobj<double>::wfobj(double *);
template wfobj<std::complex<float>>::wfobj(std::complex<float> *);
template wfobj<std::complex<double>>::wfobj(std::complex<double> *);
template wfobj<float>::~wfobj(void);
template wfobj<double>::~wfobj(void);
template wfobj<std::complex<float>>::~wfobj(void);
template wfobj<std::complex<double>>::~wfobj(void);

template spinobj<float>::spinobj(void);
template spinobj<double>::spinobj(void);
template spinobj<float>::spinobj(float *);
template spinobj<double>::spinobj(double *);
template spinobj<std::complex<float>>::spinobj(void);
template spinobj<std::complex<double>>::spinobj(void);
template spinobj<std::complex<float>>::spinobj(std::complex<float> *);
template spinobj<std::complex<double>>::spinobj(std::complex<double> *);
template spinobj<float>::~spinobj(void);
template spinobj<double>::~spinobj(void);
template spinobj<std::complex<float>>::~spinobj(void);
template spinobj<std::complex<double>>::~spinobj(void);

template void spinobj<double>::get_oppo(void);
//template void spinobj<std::complex<float>>::get_oppo(void);
//template void spinobj<std::complex<double>>::get_oppo(void);
//template void spinobj<float>::get_oppo(void);
