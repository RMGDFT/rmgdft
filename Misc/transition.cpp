// Some transitional routines


#include "rmgtypes.h"
#include "mpi.h"
#include "transition.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Subdiag.h"

BaseGrid *Rmg_G;
TradeImages *Rmg_T;
Lattice Rmg_L;

extern "C"
MPI_Comm transition_get_grid_comm(void)
{
    return pct.grid_comm;
}
extern "C" int transition_get_gridpe(void)
{
    return pct.gridpe;
}

