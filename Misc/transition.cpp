// Some transitional routines


#include "rmgtypes.h"
#include "mpi.h"
#include "transition.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "typedefs.h"

extern "C"
MPI_Comm transition_get_grid_comm(void)
{
    return pct.grid_comm;
}
