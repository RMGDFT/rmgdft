/*
 *
 * Copyright (c) 2015, Emil Briggs
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


#ifndef RMG_Pw_H
#define RMG_Pw_H 1


// This class is used to handle the integration of plane wave functionality into RMG

#include "BaseGrid.h"
#include "Lattice.h"
#include "fft3d.h"


typedef struct {
  double a[3];
} gvector;

class Pw {

private:

    // BaseGrid class
    BaseGrid *Grid;

    // Lattice object
    Lattice *L;

    // Grid spacings
    double hxgrid;
    double hygrid;
    double hzgrid;

    // true if gamma only
    bool is_gamma;

public:
    Pw (BaseGrid &G, Lattice &L, int ratio, bool gamma_flag, MPI_Comm comm);
    void index_to_gvector(int *index, double *gvector);

    ~Pw(void);

    // Real space basis on this node and globally
    int pbasis;
    int global_basis;

    // fft plans to use with this PW object
    struct fft_plan_3d *fft_forward_plan;
    struct fft_plan_3d *fft_backward_plan;

    // Real space grid dimensions on this node
    int dimx;
    int dimy;
    int dimz;

    // Global dimensions of the real space grid
    int global_dimx;
    int global_dimy;
    int global_dimz;

    // Plane wave cutoff and max
    double gcut;
    double gmax;

    int ng;
    gvector *g;
    double *gmags;
    double *gmask;

};


#endif
