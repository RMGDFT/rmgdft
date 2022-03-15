// Some transitional routines


#include "mpi.h"
#include "transition.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Subdiag.h"
#include "MpiQueue.h"

LaplacianCoeff *LC;
LaplacianCoeff *LC_6;
LaplacianCoeff *HLC;
BaseGrid *Rmg_G;
BaseGrid *Rmg_halfgrid;
TradeImages *Rmg_T;
Lattice Rmg_L;
MpiQueue *Rmg_Q;
PulayMixing *Pulay_rho;
PulayMixing *Pulay_ldau;
PulayMixing *Pulay_orbital;
