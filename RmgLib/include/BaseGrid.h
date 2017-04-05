/*
 *
 * Copyright (c) 2014, Emil Briggs
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

#ifndef RMG_BaseGrid_H
#define RMG_BaseGrid_H 1

#include "rmg_error.h"


/* Neighbor list indices */
#define NB_N 0
#define NB_S 1
#define NB_E 2
#define NB_W 3
#define NB_U 4
#define NB_D 5
#define NB_SELF 7

#ifdef __cplusplus

#include <iostream>
#include <cstdio>

/// Controls grid and nodes
class BaseGrid {

protected:
    /* Global coarse grid dimensions */
    int NX_GRID;
    int NY_GRID;
    int NZ_GRID;

    /* Total PES */
    int NPES;

    /* Node (PE) dimensions */
    int PE_X;
    int PE_Y;
    int PE_Z;

    /* Grid sizes on each PE */
    int PX0_GRID;
    int PY0_GRID;
    int PZ0_GRID;

    /* Basis size on each PE */
    int P0_BASIS;

    /* MPI specific info */
    int rank;
    int neighbors[6];
    bool face_pe;
    bool edge_pe;
    bool corner_pe;
  
private:

    /* Grid offsets on each PE */
    int PX_OFFSET;
    int PY_OFFSET;
    int PZ_OFFSET;

public:

    BaseGrid(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int rank, int default_FG_RATIO);

    /* Fine grid/coarse default ratio */
    int default_FG_RATIO;

    /* Function prototypes */
    void set_rank(int new_rank, MPI_Comm comm);
    void find_node_sizes(int rank, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize);
    void find_node_offsets(int rank, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset);

    int get_default_FG_RATIO(void);

    int get_NPES(void);
    int get_PE_X(void);
    int get_PE_Y(void);
    int get_PE_Z(void);

    int get_NX_GRID(int density);
    int get_NY_GRID(int density);
    int get_NZ_GRID(int density);

    double get_hxgrid(int density);
    double get_hygrid(int density);
    double get_hzgrid(int density);

    int get_PX0_GRID(int density);
    int get_PY0_GRID(int density);
    int get_PZ0_GRID(int density);

    int get_PX_OFFSET(int density);
    int get_PY_OFFSET(int density);
    int get_PZ_OFFSET(int density);

    int get_P0_BASIS(int density);
    int get_GLOBAL_BASIS(int density);

    void pe2xyz(int pe, int *x, int *y, int *z);
    int xyz2pe(int x, int y, int z);

    bool is_face_pe(void);
    bool is_edge_pe(void);
    bool is_corner_pe(void);

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *get_neighbors(void);

    // Returns the rank of this process
    int get_rank(void);
    MPI_Comm comm;

};

#endif
#endif
