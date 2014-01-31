

#include "const.h"
#include "grid.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "mg.h"
#include <math.h>

int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxsize, int *pysize, int *pzsize)
{

    int ii, jj, kk;
    int idx, ix, iy, iz, mfac;

    pe2xyz (gridpe, &ii, &jj, &kk);

    mfac = nxgrid / NX_GRID;
    *pxsize = mfac * (NX_GRID/PE_X);
    ix = NX_GRID % PE_X;
    if(ii < ix) *pxsize += mfac;

    mfac = nygrid / NY_GRID;
    *pysize = mfac * (NY_GRID/PE_Y);
    iy = NY_GRID % PE_Y;
    if(jj < iy) *pysize += mfac;

    mfac = nzgrid / NZ_GRID;
    *pzsize = mfac * (NZ_GRID/PE_Z);
    iz = NZ_GRID % PE_Z;
    if(kk < iz) *pzsize += mfac;
}


