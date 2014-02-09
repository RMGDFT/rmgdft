
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"
#include <math.h>

// Finds a 1-D node coordinate that a given coarse grid point
// falls into.
// igrid = 1-D index of global coarse grid point
// tgrid = dimension of global coarse grid point in same direction
// pgrid = node (PE) dimensions in same direction
void find_grid_owner(int igridx, int igridy, int igridz, int nxgrid, int nygrid, int nzgrid, int *xnode, int *ynode, int *znode)
{
    int mingrid, maxgrid, remainder, mfac;
 
    mfac = nxgrid / get_NX_GRID();
    mingrid = mfac * (get_NX_GRID() / get_PE_X());
    maxgrid = mingrid;
    remainder = get_NX_GRID() % get_PE_X();
    if(remainder) maxgrid += mfac;
    if(maxgrid == mingrid) {
        *xnode = igridx / maxgrid;
    }
    else {
        if(igridx < (remainder * maxgrid)) {
            *xnode = igridx / maxgrid;
        }
        else {
            *xnode = remainder;
            igridx -= remainder * maxgrid;
            *xnode += igridx /mingrid;
        }
    }

    mfac = nygrid / get_NY_GRID();
    mingrid = mfac * (get_NY_GRID() / get_PE_Y());
    maxgrid = mingrid;
    remainder = get_NY_GRID() % get_PE_Y();
    if(remainder) maxgrid += mfac;
    if(maxgrid == mingrid) {
        *ynode = igridy / maxgrid;
    }
    else {
        if(igridy < remainder * maxgrid) {
            *ynode = igridy / maxgrid;
        }
        else {
            *ynode = remainder;
            igridy -= remainder * maxgrid;
            *ynode += igridy /mingrid;
        }
    }

    mfac = nzgrid / get_NZ_GRID();
    mingrid = mfac * (get_NZ_GRID() / get_PE_Z());
    maxgrid = mingrid;
    remainder = get_NZ_GRID() % get_PE_Z();
    if(remainder) maxgrid += mfac;
    if(maxgrid == mingrid) {
        *znode = igridz / maxgrid;
    }
    else {
        if(igridz < remainder * maxgrid) {
            *znode = igridz / maxgrid;
        }
        else {
            *znode = remainder;
            igridz -= remainder * maxgrid;
            *znode += igridz /mingrid;
        }
    }

}
