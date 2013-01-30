

#include "main.h"
#include <math.h>

int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset)
{

    int ii, jj, kk;
    int idx, ix, iy, iz, ioffset;

    pe2xyz (gridpe, &ii, &jj, &kk);

    // Now compute the global grid offset of the first point of the node grid
    *pxoffset = ii*(nxgrid/PE_X);
    ix = nxgrid % PE_X;
    ioffset = 0;
    for(idx = 0;idx < ii;idx++) {
        if(ix && (idx < ix)) ioffset++;
    }
    *pxoffset += ioffset;

    *pyoffset = jj*(nygrid/PE_Y);
    iy = nygrid % PE_Y;
    ioffset = 0;
    for(idx = 0;idx < jj;idx++) {
        if(iy && (idx < iy)) ioffset++;
    }
    *pyoffset += ioffset;

    *pzoffset = kk*(nzgrid/PE_Z);
    iz = nzgrid % PE_Z;
    ioffset = 0;
    for(idx = 0;idx < kk;idx++) {
        if(iz && (idx < iz)) ioffset++;
    }
    *pzoffset += ioffset;

}


