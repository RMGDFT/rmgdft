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
     Same as above except existing storage at *data is used and is not deallocated
     when the object goes out of scope. Responsibility for that is left with the
     upper level routines.

     GridObject should never be created directly though. Instead use the derived
     classes which set the size and other characteristics correctly for the
     type of object in question.

     For potentials like vh, vnuc on the fine grid
         fgobj() V;
         fgobj() V(ptr);

     For spin dependent objects on the fine grid. The number of components for these
     varies depending on the type of calculation
     1 for non spin
     2 for spin
     4 for spin orbit
         spinobj() rho;
         spinobj() rho(ptr);
         spinobj() vxc;
         spinobj() vxc(ptr);

     Wavefunctions are defined on the coarse grid with 1 or 2 components.
     1 for non-spin and spin polarized
     2 for spin orbit
         wfobj() x;
         x.up[] or X.dw[]

     
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

    T& operator [](int idx) {
        return data_[idx];
    }

    T operator [](int idx) const {
        return data_[idx];
    }

private:
   bool owns_allocation = true;

protected:
   int dimx_;
   int dimy_;
   int dimz_;
   int pbasis;
   int factor = 1;
   T *data_;

   void allocate(int components)
   {
       factor = components;
       data_ = new T[factor*pbasis]();
   }
   void allocate(int components, T *ptr)
   {
       factor = components;
       data_ = ptr;
       owns_allocation = false;
   }
   void increment( const GridObject& c );
   void decrement( const GridObject& c );
   void multiply( const T& b );

};


template <typename T> class fgobj : public GridObject<T>
{
public:
    fgobj(void);
    fgobj(T *data_ptr);
    ~fgobj(void);
};

template <typename T> class spinobj : public GridObject<T>
{
public:
    spinobj(void);
    spinobj(T *data_ptr);
    ~spinobj(void);
    std::span<T> up;
    std::span<T> dw;
    std::span<T> c0;
    std::span<T> cx;
    std::span<T> cy;
    std::span<T> cz;

};

// For spin-orbit there are two components for each wavefunction
// references as up or dw.
template <typename T> class wfobj : public GridObject<T>
{
public:
    wfobj(void);
    wfobj(T *data_ptr);
    ~wfobj(void);
    std::span<T> up;
    std::span<T> dw;
};


#endif
