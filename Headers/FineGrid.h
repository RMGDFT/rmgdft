#ifndef RMG_FineGrid_H
#define RMG_FineGrid_H 1

#include "BaseGrid.h"

// Uncomment when generating doxygen docs
//#define __cplusplus
#ifdef __cplusplus

/// Controls grid and nodes
class FineGrid : public BaseGrid {

private:
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
