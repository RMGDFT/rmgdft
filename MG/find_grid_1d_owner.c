
#include "main.h"
#include <math.h>

// Finds a 1-D node coordinate that a given coarse grid point
// falls into.
// igrid = 1-D index of global coarse grid point
// tgrid = dimension of global coarse grid point in same direction
// pgrid = node (PE) dimensions in same direction
int find_grid_1d_owner(int igrid, int tgrid, int pgrid)
{
    int mingrid, maxgrid, remainder;
 
    mingrid = tgrid / pgrid; 
    maxgrid = mingrid;
    remainder = tgrid % pgrid;
    if(remainder) maxgrid++;

    if(maxgrid == mingrid) return igrid / mingrid;

    if(igrid < (remainder * maxgrid)) return igrid / maxgrid;
    return (remainder + igrid / mingrid);
}
