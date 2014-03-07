#ifndef RMG_BaseGrid_H
#define RMG_BaseGrid_H 1

#include "rmgtypes.h"
#include "rmg_error.h"


/* Neighbor list indices */
#define NB_N 0
#define NB_S 1
#define NB_E 2
#define NB_W 3
#define NB_U 4
#define NB_D 5
#define NB_SELF 7

// Uncomment when generating doxygen docs
//#define __cplusplus
#ifdef __cplusplus

#include <iostream>
#include <cstdio>

/// Controls grid and nodes
class BaseGrid {

protected:
    /* Global coarse grid dimensions */
    static int NX_GRID;
    static int NY_GRID;
    static int NZ_GRID;

    /* Node (PE) dimensions */
    static int PE_X;
    static int PE_Y;
    static int PE_Z;

    /* Grid sizes on each PE */
    static int PX0_GRID;
    static int PY0_GRID;
    static int PZ0_GRID;

    /* Basis size on each PE */
    static int P0_BASIS;

    /* Fine grid basis size on each PE */
    static int FP0_BASIS;

    /* MPI specific info */
    static int gridpe;
    static int neighbors[6];

private:

    /* Global fine grid dimensions */
    static int FNX_GRID;
    static int FNY_GRID;
    static int FNZ_GRID;

    /* Grid offsets on each PE */
    static int PX_OFFSET;
    static int PY_OFFSET;
    static int PZ_OFFSET;

    /* Fine grid sizes on each PE */
    static int FPX0_GRID;
    static int FPY0_GRID;
    static int FPZ0_GRID;

    /* Fine Grid offsets on each PE */
    static int FPX_OFFSET;
    static int FPY_OFFSET;
    static int FPZ_OFFSET;

    /* Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. A value larger than 1.05 can lead to convergence problems. */
    static rmg_double_t anisotropy;

    /* Initialiazation flags */
    static int neighbor_first;
    static int grid_first;
    static int anisotropy_first;

public:

    /* Fine grid/coarse grid ratio */
    static int FG_NX;
    static int FG_NY;
    static int FG_NZ;

    /* Function prototypes */
    void set_grids(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int FG_NX, int FG_NY, int FG_NZ);
    void set_nodes(int newgridpe, int ii, int jj, int kk);
    void set_neighbors(int *list);
    int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize);
    int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset);

    int get_PE_X(void);
    int get_PE_Y(void);
    int get_PE_Z(void);

    int get_NX_GRID(void);
    int get_NY_GRID(void);
    int get_NZ_GRID(void);

    int get_FNX_GRID(void);
    int get_FNY_GRID(void);
    int get_FNZ_GRID(void);

    int get_PX0_GRID(void);
    int get_PY0_GRID(void);
    int get_PZ0_GRID(void);

    int get_PX_OFFSET(void);
    int get_PY_OFFSET(void);
    int get_PZ_OFFSET(void);

    int get_FPX_OFFSET(void);
    int get_FPY_OFFSET(void);
    int get_FPZ_OFFSET(void);

    int get_P0_BASIS(void);
    int get_FP0_BASIS(void);

    int get_FPX0_GRID(void);
    int get_FPY0_GRID(void);
    int get_FPZ0_GRID(void);

    void set_anisotropy(rmg_double_t a);
    rmg_double_t get_anisotropy(void);
    void pe2xyz(int pe, int *x, int *y, int *z);

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *get_neighbors(void);

    // Returns the rank of this process
    int get_gridpe(void);

};

#endif
#endif
