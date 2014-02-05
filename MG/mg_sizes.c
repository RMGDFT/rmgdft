
#include "const.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "mg.h"
#include "main.h"
#include <float.h>
#include <math.h>


/* Compute 1-D grid sizes for the next multigrid level 

Inputs:
curdim        = current size of this grid on this node
global_dim    = global grid dimension
global_offset = offset of edge of this node grid on the global grid
global_pdim   = dimension of this node grid
bctype        = boundary condition


Outputs:
*roffset      = pointer to grid offset (always 0 or 1)

Return value  = size of next grid level


*/


int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
{
    int skip, new_dim, istart, istop;

    // Default offset is 0
    *roffset = 0;

    if(bctype == PERIODIC) {

        skip = (2 << curlevel);
        // First check if we have too many multigrid levels. For periodic boundary
        // conditions the next level of the global grid must be divisible by 2
        if ((global_dim % skip) != 0) {
            error_handler ("Too many multigrid levels specified.");
        }

        // Require at least one point in the level
        new_dim = global_pdim / skip;
        if(!new_dim) {
            error_handler ("Too many multigrid levels specified.");
        }

        // evenly divisible then we are done
        if(!(global_pdim % skip)) return new_dim;

        // Check if first point is included and if not subtract
        istart = skip - global_offset % skip;
        istop = (global_offset + global_pdim - 1) % skip;
        if((istart == skip) || (istop == skip)) new_dim++;
        
        // Perform offset check
        if((istart == skip) || (istart == 0)) {
            return new_dim;
        }
        *roffset = 1;
        
        return new_dim;

    }

    rmg_error_handler("Boundary condition not programmed."); 

}

