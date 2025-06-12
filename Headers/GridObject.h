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

template <typename T> class fgobj;
template <typename T> class wfobj;
template <typename T> class spinobj;

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
    return a;
  }
  friend GridObject operator*(GridObject a, const T &b) {
    a *= b;
    return a;
  }
  friend GridObject& operator*=(GridObject& a, const T &b) {
    a.multiply(b);
    return a;
  }

public:
    GridObject(int density);
    GridObject(int density, T *data_ptr);
    GridObject(GridObject& t) // copy constructor
    {
        this->dimx_  = t.dimx;
        this->dimy_  = t.dimy;
        this->dimz_  = t.dimz;
        this->incx_  = t.incx;
        this->incy_  = t.incy;
        this->incz_  = t.incz;
        this->offsetx_  = t.offsetx;
        this->offsety_  = t.offsety;
        this->offsetz_  = t.offsetz;
        this->pbasis_   = t.pbasis;
        this->factor    = t.factor;
        this->allocate(this->factor);
        this->owns_allocation = true;
        std::copy(t.data_, t.data_ + this->factor*this->pbasis_, this->data_);
    }
    ~GridObject(void);

    void set(T x)
    {
        for(int i=0;i < this->factor*this->pbasis_;i++) this->data_[i] = x;
    }

    // Dimensions and offsets on each MPI task
    // These are public and read only external to the
    // class but reference the internally writeable data
    const int& dimx = dimx_;
    const int& dimy = dimy_;
    const int& dimz = dimz_;
    const int& incx = incx_;
    const int& incy = incy_;
    const int& incz = incz_;
    const int& offsetx = offsetx_;
    const int& offsety = offsety_;
    const int& offsetz = offsetz_;
    const int& pbasis = pbasis_;
    const double& vel = vel_;

    const int size() const { return pbasis; }
    T* data() { return data_; }

    T& operator [](int idx) {
        return data_[idx];
    }

    T operator [](int idx) const {
        return data_[idx];
    }

    T& operator ()(int i, int j, int k) {
        return data_[i*incx_ + j*incy_ + k];
    }

    T operator ()(int i, int j, int k) const {
        return data_[i*incx_ + j*incy_ + k];
    }


private:
   bool owns_allocation = true;

protected:
   int dimx_;
   int dimy_;
   int dimz_;
   int incx_;
   int incy_;
   int incz_;
   int offsetx_;
   int offsety_;
   int offsetz_;
   int pbasis_;
   int factor = 1;
   double vel_;
   T *data_;

   void allocate(int components)
   {
       factor = components;
       data_ = new T[factor*pbasis_]();
   }
   void allocate(int components, T *ptr)
   {
       factor = components;
       data_ = ptr;
       owns_allocation = false;
   }
   void increment( const GridObject& c );
   void increment( const fgobj<T>& c );
   void increment( const wfobj<T>& c );
   void increment( const spinobj<T>& c );
   void decrement( const GridObject& c );
   void decrement( const fgobj<T>& c );
   void decrement( const wfobj<T>& c );
   void decrement( const spinobj<T>& c );
   void multiply( const T& b );

   GridObject& operator=(GridObject const &rhs)
   {
       if (this != &rhs)
       {
           std::copy(rhs.data_, rhs.data_ + this->factor*this->pbasis_, this->data_);
       }
       return *this;
    };

};


template <typename T> class fgobj : public GridObject<T>
{

  friend fgobj operator+(fgobj a, const fgobj& b) {
    a += b;
    return std::move(a);
  }
  friend fgobj& operator+=(fgobj& a, const fgobj& b) {
    a.increment(b);
    return a;
  }
  friend fgobj operator-(fgobj a, const fgobj& b) {
    a -= b;
    return std::move(a);
  }
  friend fgobj& operator-=(fgobj& a, const fgobj& b) {
    a.decrement(b);
    return a;
  }
  friend fgobj operator*(const T &b, fgobj a) {
    a.multiply(b);
    return a;
  }
  friend fgobj operator*(fgobj a, const T &b) {
    a *= b;
    return a;
  }
  friend fgobj& operator*=(fgobj& a, const T &b) {
    a.multiply(b);
    return a;
  }

public:
    fgobj(void);
    fgobj(T *data_ptr);
    ~fgobj(void);

protected:
//   void increment( const fgobj& c );
//   void decrement( const fgobj& c );

};

template <typename T> class spinobj : public GridObject<T>
{
  friend spinobj operator+(spinobj a, const spinobj& b) {
    a += b;
    return a;
  }
  friend spinobj& operator+=(spinobj& a, const spinobj& b) {
    a.increment(b);
    return a;
  }
  friend spinobj operator-(spinobj a, const spinobj& b) {
    a -= b;
    return a;
  }
  friend spinobj& operator-=(spinobj& a, const spinobj& b) {
    a.decrement(b);
    return a;
  }
  friend spinobj operator*(const T &b, spinobj a) {
    a.multiply(b);
    return a;
  }
  friend spinobj operator*(spinobj a, const T &b) {
    a *= b;
    return a;
  }
  friend spinobj& operator*=(spinobj& a, const T &b) {
    a.multiply(b);
    return a;
  }
public:
    spinobj(void);
    spinobj(T *data_ptr);
    ~spinobj(void);
    void get_oppo(void);
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
  friend wfobj operator+(wfobj a, const wfobj& b) {
    a += b;
    return std::move(a);
  }
  friend wfobj& operator+=(wfobj& a, const wfobj& b) {
    a.increment(b);
    return a;
  }
  friend wfobj operator-(wfobj a, const wfobj& b) {
    a -= b;
    return std::move(a);
  }
  friend wfobj& operator-=(wfobj& a, const wfobj& b) {
    a.decrement(b);
    return a;
  }
  friend wfobj operator*(const T &b, wfobj a) {
    a.multiply(b);
    return a;
  }
  friend wfobj operator*(wfobj a, const T &b) {
    a *= b;
    return a;
  }
  friend wfobj& operator*=(wfobj& a, const T &b) {
    a.multiply(b);
    return a;
  }
public:
    wfobj(void);
    wfobj(T *data_ptr);
    ~wfobj(void);
    std::span<T> up;
    std::span<T> dw;
};


#endif
