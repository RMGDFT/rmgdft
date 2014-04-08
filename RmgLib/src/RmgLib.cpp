
#include "RmgLib.h"

/// Global coarse grid X dimension
int RmgLib::NX_GRID;
// Global coarse grid Y dimension
int RmgLib::NY_GRID;
/// Global coarse grid Z dimension
int RmgLib::NZ_GRID;

/* Node (PE) dimensions */
int RmgLib::PE_X;
int RmgLib::PE_Y;
int RmgLib::PE_Z;

int RmgLib::default_FG_RATIO;

// initialization flag
bool RmgLib::initialized = false;

/// Used to set up global coarse grid dimensions, MPI node dimensions and the ratio of the fine grid to the coarse in each coordinate dimension.
/// @param newNX_GRID New global coarse grid X dimension
/// @param newNY_GRID New global coarse grid Y dimension
/// @param newNZ_GRID New global coarse grid Z dimension
/// @param newPE_X New MPI grid X dimension
/// @param newPE_Y New MPI grid Y dimension
/// @param newPE_Z New MPI grid Z dimension
/// @param newFG_RATIO New ratio of fine grid to coarse
void RmgLib::set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_RATIO)
{

    RmgLib::NX_GRID = newNX_GRID;
    RmgLib::NY_GRID = newNY_GRID;
    RmgLib::NZ_GRID = newNZ_GRID;

    RmgLib::PE_X = newPE_X;
    RmgLib::PE_Y = newPE_Y;
    RmgLib::PE_Z = newPE_Z;

    RmgLib::default_FG_RATIO = newFG_RATIO;

    RmgLib::initialized = true;
}
