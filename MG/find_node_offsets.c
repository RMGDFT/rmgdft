

#include "main.h"
#include <math.h>

int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset)
{

    int ii, jj, kk;
    int idx, ix, iy, iz, ioffset, mfac;

    pe2xyz (gridpe, &ii, &jj, &kk);

    // Now compute the global grid offset of the first point of the node grid
    mfac = nxgrid / NX_GRID;
    *pxoffset = mfac * ii*(NX_GRID/PE_X);
    ix = NX_GRID % PE_X;
    ioffset = 0;
    for(idx = 1;idx <= ii;idx++) {
        if(idx <= ix) ioffset++;
    }
    ioffset *= mfac;
    *pxoffset += ioffset;


    mfac = nygrid / NY_GRID;
    *pyoffset = mfac * jj*(NY_GRID/PE_Y);
    iy = NY_GRID % PE_Y;
    ioffset = 0;
    for(idx = 1;idx <= jj;idx++) {
        if(idx <= iy) ioffset++;
    }
    ioffset *= mfac;
    *pyoffset += ioffset;

    mfac = nzgrid / NZ_GRID;
    *pzoffset = mfac * kk*(NZ_GRID/PE_Z);
    iz = NZ_GRID % PE_Z;
    ioffset = 0;
    for(idx = 1;idx <= kk;idx++) {
        if(idx <= iz) ioffset++;
    }
    ioffset *= mfac;
    *pzoffset += ioffset;

}


