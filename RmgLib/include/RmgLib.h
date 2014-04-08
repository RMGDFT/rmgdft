#ifndef RMG_RmgLib_H
#define RMG_RmgLib_H 1

/// Controls grid and nodes
class RmgLib {

protected:

    /* Global coarse grid dimensions */
    static int NX_GRID;
    static int NY_GRID;
    static int NZ_GRID;

    /* Node (PE) dimensions */
    static int PE_X;
    static int PE_Y;
    static int PE_Z;

    static bool initialized;

    static int default_FG_RATIO;

public:
    void set_grids(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int default_FG_RATIO);


};

#endif

