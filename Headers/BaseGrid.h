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

    /* MPI specific info */
    static int gridpe;
    static int neighbors[6];

private:

    /* Grid offsets on each PE */
    static int PX_OFFSET;
    static int PY_OFFSET;
    static int PZ_OFFSET;

    /* Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. A value larger than 1.05 can lead to convergence problems. */
    static double anisotropy;

    /* Initialiazation flags */
    static int neighbor_first;
    static int grid_first;
    static int anisotropy_first;

public:

    /* Fine grid/coarse default ratio */
    static int default_FG_RATIO;

    /* Function prototypes */
    void set_grids(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int default_FG_RATIO);
    void set_nodes(int newgridpe);
    void find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize);
    void find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset);

    int get_default_FG_RATIO(void);

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

    void set_anisotropy(double a);
    double get_anisotropy(void);
    void pe2xyz(int pe, int *x, int *y, int *z);
    int xyz2pe(int x, int y, int z);

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *get_neighbors(void);

    // Returns the rank of this process
    int get_gridpe(void);

};

#endif
#endif
