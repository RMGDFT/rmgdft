

#include "main.h"
#include "common_prototypes.h"
#include "grid.h"
#include <math.h>

int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset)
{

    int ii, jj, kk;
    int idx, ix, iy, iz, ioffset, mfac;

    pe2xyz (gridpe, &ii, &jj, &kk);

    // Now compute the global grid offset of the first point of the node grid
    mfac = nxgrid / get_NX_GRID();
    *pxoffset = mfac * ii*(get_NX_GRID()/get_PE_X());
    ix = get_NX_GRID() % get_PE_X();
    ioffset = 0;
    for(idx = 1;idx <= ii;idx++) {
        if(idx <= ix) ioffset++;
    }
    ioffset *= mfac;
    *pxoffset = *pxoffset + ioffset;


    mfac = nygrid / get_NY_GRID();
    *pyoffset = mfac * jj*(get_NY_GRID()/get_PE_Y());
    iy = get_NY_GRID() % get_PE_Y();
    ioffset = 0;
    for(idx = 1;idx <= jj;idx++) {
        if(idx <= iy) ioffset++;
    }
    ioffset *= mfac;
    *pyoffset += ioffset;

    mfac = nzgrid / get_NZ_GRID();
    *pzoffset = mfac * kk*(get_NZ_GRID()/get_PE_Z());
    iz = get_NZ_GRID() % get_PE_Z();
    ioffset = 0;
    for(idx = 1;idx <= kk;idx++) {
        if(idx <= iz) ioffset++;
    }
    ioffset *= mfac;
    *pzoffset += ioffset;

}


