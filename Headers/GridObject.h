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

*/

#ifndef RMG_GridObject_H
#define RMG_GridObject_H 1


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
    GridObject(int density, T *p);
    ~GridObject(void);
    int dimx() const { return dimx_; }
    int dimy() const { return dimy_; }
    int dimz() const { return dimz_; }
    T* data() { return data_.data(); }

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
   std::vector<T> data_;
   bool owns_allocation = true;

protected:
  void increment( const GridObject& c );
  void decrement( const GridObject& c );
  void multiply( const T& b );

};

template<typename T>
GridObject<T>::GridObject(int density)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    pbasis = dimx_ * dimy_ * dimz_;
    data_.resize(pbasis);
}

template<typename T>
GridObject<T>::GridObject(int density, T *data)
{
    dimx_ = Rmg_G->get_PX0_GRID(density);
    dimy_ = Rmg_G->get_PY0_GRID(density);
    dimz_ = Rmg_G->get_PZ0_GRID(density);
    pbasis = dimx_ * dimy_ * dimz_;
    data_ = data;
    owns_allocation = false;
}

template<typename T>
GridObject<T>::~GridObject(void)
{
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
#endif
