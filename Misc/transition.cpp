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

// Pointer to Kpoint class array
/* Hartree potential */
extern double *vh;
/* Nuclear local potential */
extern double *vnuc;
/* Exchange-correlation potential */
extern double *vxc;
// Pointer to Kpoint class array
extern void **Kptr;

extern "C"
MPI_Comm transition_get_grid_comm(void)
{
    return pct.grid_comm;
}
extern "C" int transition_get_gridpe(void)
{
    return pct.gridpe;
}

extern "C" void subdiag_gamma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
     Kpoint<double>* kptr = (Kpoint<double> *)Kptr[0]; 
     Subdiag<double> (kptr, vh, vnuc, vxc, ct.subdiag_driver);
}

