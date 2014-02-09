

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

    mfac = nxgrid / get_NX_GRID();
    *pxsize = mfac * (get_NX_GRID()/get_PE_X());
    ix = get_NX_GRID() % get_PE_X();
    if(ii < ix) *pxsize += mfac;

    mfac = nygrid / get_NY_GRID();
    *pysize = mfac * (get_NY_GRID()/get_PE_Y());
    iy = get_NY_GRID() % get_PE_Y();
    if(jj < iy) *pysize += mfac;

    mfac = nzgrid / get_NZ_GRID();
    *pzsize = mfac * (get_NZ_GRID()/get_PE_Z());
    iz = get_NZ_GRID() % get_PE_Z();
    if(kk < iz) *pzsize += mfac;
}


