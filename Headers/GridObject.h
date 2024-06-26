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


/* This class implements GridObjects for RMG. A GridObject encapsulates something stored on
   a real space domain such as wavefunctions, charge density etc. The object knows it's size
   and topology so working with it will be simpler and cleaner than working with the lower
   level data. The class includes operator overloads that make it easier to perform standard
   operations.


   Constructor examples:
     GridObject<T> V(density);
     Creates a grid object V on the grid specified by density
     (1= wavefunction grid, >1 = charge density grid). Storage
     for the object is allocated internally and destroyed when
     the object is deleted or goes out of scope.


     GridObject<T> V(density, T *data);
     Same as above except existing storage at *data is used. 

     One can also create GridObjects specific to the fine or the coarse Grid.
     FineGridObject() V;
     FineGridObject() V(ptr);
     WfGridObject() V;
     WfGridObject() V(ptr);

     For spin orbit (ct.noncoll_factor=2) there are two components to a WfGridObject
     referenced as up or dw.
     WfGridObject() X;
     x.up[] or X.dw[]

     There is also a SpinFineGridObject which is identical to FineGridObject for
     non spin polarized calculations but has two components when ct.nspin = 2.

*/

#ifndef RMG_GridObject_H
#define RMG_GridObject_H 1

#include <complex>
#include <span>


template <typename T> class GridObject {

  friend GridObject operator+(GridObject a, const GridObject& b) {
    a += b;
    return std::move(a);
  }
  friend GridObject& operator+=(GridObject& a, const GridObject& b) {
    a.increment(b);
    return a;
  }
  friend GridObject operator-(GridObject a, const GridObject& b) {
    a -= b;
    return std::move(a);
  }
  friend GridObject& operator-=(GridObject& a, const GridObject& b) {
    a.decrement(b);
    return a;
  }
  friend GridObject operator*(const T &b, GridObject a) {
    a.multiply(b);
    return std::move(a);
  }
  friend GridObject operator*(GridObject a, const T &b) {
    a.multiply(b);
    return std::move(a);
  }
  friend GridObject& operator*=(GridObject& a, const T &b) {
    a.multiply(b);
    return a;
  }

public:
    GridObject(int density);
    GridObject(int density, T *data_ptr);
    ~GridObject(void);
    int dimx() const { return dimx_; }
    int dimy() const { return dimy_; }
    int dimz() const { return dimz_; }
    const int size() const { return pbasis; }
    T* data() { return data_; }
    std::span<T> up;
    std::span<T> dw;

    T& operator [](int idx) {
        return data_[idx];
    }

    T operator [](int idx) const {
        return data_[idx];
    }

private:
   int dimx_;
   int dimy_;
   int dimz_;
   int pbasis;
   T *data_;
   bool owns_allocation = true;

protected:
   void allocate(int components)
   {
       factor = components;
       data_ = new T[factor*pbasis]();
       up = std::span(data_, pbasis);
       dw = std::span(data_ + (factor-1)*pbasis, pbasis);
   }
   void allocate(int components, T *ptr)
   {
       factor = components;
       data_ = ptr;
       up = std::span(data_, pbasis);
       dw = std::span(data_ + (factor-1)*pbasis, pbasis);
       owns_allocation = false;
   }

   int factor = 1;
   void increment( const GridObject& c );
   void decrement( const GridObject& c );
   void multiply( const T& b );

};


template <typename T> class FineGridObject : public GridObject<T>
{
public:
    FineGridObject(void);
    FineGridObject(T *data_ptr);
    ~FineGridObject(void);
#if 0
    FineGridObject(void) : GridObject<T>(Rmg_G->default_FG_RATIO)
    {
        this->allocate(1);
    }
    FineGridObject(T *data_ptr) : GridObject<T>(Rmg_G->default_FG_RATIO, data_ptr)
    {
        this->allocate(1, data_ptr);
    }
    ~FineGridObject(void)
    {
    }
#endif
};

template <typename T> class SpinFineGridObject : public GridObject<T>
{
public:
    SpinFineGridObject(void);
    SpinFineGridObject(T *data_ptr);
    ~SpinFineGridObject(void);
#if 0
    SpinFineGridObject(void) : GridObject<T>(Rmg_G->default_FG_RATIO)
    {
        this->allocate(ct.nspin);
    }
    SpinFineGridObject(T *data_ptr) : GridObject<T>(Rmg_G->default_FG_RATIO, data_ptr)
    {
        this->allocate(ct.nspin, data_ptr);
    }
    ~SpinFineGridObject(void)
    {
    }
#endif
};

// For spin-orbit there are two components for each wavefunction
// references as up or dw.
template <typename T> class WfGridObject : public GridObject<T>
{
public:
    WfGridObject(void) : GridObject<T>(1)
    {
        this->allocate(ct.noncoll_factor);
    }
    WfGridObject(T *data_ptr) : GridObject<T>(1, data_ptr)
    {
        this->allocate(ct.noncoll_factor, data_ptr);
    }
    ~WfGridObject(void)
    {
    }
};


#endif
