/**
 * @file
 *
 *
 * @section DESCRIPTION
 * Used to access grid and node related data.
 */

#include "FineGrid.h"


using namespace std;

int GLOBAL_GRIDX;
int GLOBAL_GRIDY;
int GLOBAL_GRIDZ;
int PE_GRIDX;
int PE_GRIDY;
int PE_GRIDZ;
int PE_OFFSETX;
int PE_OFFSETY;
int PE_OFFSETZ;
int GLOBAL_BASIS;
int PE_BASIS;

    FineGrid::FineGrid(int density)
    {
        GLOBAL_GRIDX = density * BaseGrid::NX_GRID; 
        GLOBAL_GRIDY = density * BaseGrid::NY_GRID; 
        GLOBAL_GRIDZ = density * BaseGrid::NZ_GRID; 
        PE_GRIDX = density * BaseGrid::PX0_GRID; 
        PE_GRIDY = density * BaseGrid::PY0_GRID; 
        PE_GRIDZ = density * BaseGrid::PZ0_GRID; 
        PE_BASIS = density * density * density * BaseGrid::P0_BASIS;
        GLOBAL_BASIS = density * density * density * BaseGrid::FP0_BASIS;

        BaseGrid::find_node_offsets(BaseGrid::gridpe, GLOBAL_GRIDX, GLOBAL_GRIDY, GLOBAL_GRIDZ,
                          &PE_OFFSETX, &PE_OFFSETY, &PE_OFFSETZ);

    }

    int FineGrid::get_GLOBAL_GRIDX(void)
    {
        return GLOBAL_GRIDX;
    }                             

    int FineGrid::get_GLOBAL_GRIDY(void)
    {
        return GLOBAL_GRIDY;
    }

    int FineGrid::get_GLOBAL_GRIDZ(void)
    {
        return GLOBAL_GRIDZ;
    }

    int FineGrid::get_PE_GRIDX(void)
    {
        return PE_GRIDX;
    }                             

    int FineGrid::get_PE_GRIDY(void)
    {
        return PE_GRIDY;
    }

    int FineGrid::get_PE_GRIDZ(void)
    {
        return PE_GRIDZ;
    }

    int FineGrid::get_PE_OFFSETX(void)
    {
        return PE_OFFSETX;
    }

    int FineGrid::get_PE_OFFSETY(void)
    {
        return PE_OFFSETY;
    }

    int FineGrid::get_PE_OFFSETZ(void)
    {
        return PE_OFFSETZ;
    }

    int FineGrid::get_GLOBAL_BASIS(void)
    {
        return GLOBAL_BASIS;
    }

    int FineGrid::get_PE_BASIS(void)
    {
        return PE_BASIS;
    }
