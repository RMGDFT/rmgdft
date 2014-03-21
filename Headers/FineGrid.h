#ifndef RMG_FineGrid_H
#define RMG_FineGrid_H 1

#include "BaseGrid.h"

// Uncomment when generating doxygen docs
//#define __cplusplus
#ifdef __cplusplus

/// Controls grid and nodes
class FineGrid : public BaseGrid {

private:
    // Node (PE) dimensions are inherited from BaseGrid
    // static int PE_X;
    // static int PE_Y;
    // static int PE_Z;

    // Global coarse grid dimensions (NX_GRID, NY_GRID, NZ_GRID) are inherited from BaseGrid
    // Global dimensions for this grid
    int GLOBAL_GRIDX;
    int GLOBAL_GRIDY;
    int GLOBAL_GRIDZ;

    // Per node grid dimensions for this grid
    int PE_GRIDX;
    int PE_GRIDY;
    int PE_GRIDZ;

    // Multiplier is the integer ratio of this grid to the coarse grid
    /* Grid sizes on each PE */
    int GRID_MULTIPLIER;

    /* Grid offsets on each PE */
    int PE_OFFSETX;
    int PE_OFFSETY;
    int PE_OFFSETZ;

    // Global basis size
    int GLOBAL_BASIS;

    // Basis size on each PE 
    int PE_BASIS;

public:

    /* Function prototypes */
//    void set_grids(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int FG_NX, int FG_NY, int FG_NZ);
//    void set_nodes(int newgridpe, int ii, int jj, int kk);
//    int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize);
//    int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset);

    FineGrid(int density);

    int get_GLOBAL_GRIDX(void);
    int get_GLOBAL_GRIDY(void);
    int get_GLOBAL_GRIDZ(void);

    int get_PE_GRIDX(void);
    int get_PE_GRIDY(void);
    int get_PE_GRIDZ(void);

    int get_PE_OFFSETX(void);
    int get_PE_OFFSETY(void);
    int get_PE_OFFSETZ(void);

    int get_GLOBAL_BASIS(void);
    int get_PE_BASIS(void);

};

#endif
#endif
