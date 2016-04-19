/************************** SVN Revision Information **************************
 **    $Id: UpdatePot.cpp 3198 2015-09-02 21:27:36Z luw $    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"

#include "prototypes_on.h"

void VhDriver(double *rho_tot, double *rhoc, double *vh)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int idx, ione = 1;


    if(ct.poisson_solver == POISSON_PFFT_SOLVER)
    {
        RmgTimer *RT2 = new RmgTimer("VhPfft");
        VhPfft(rho_tot, rhoc, vh);
        delete(RT2);
    }
    else
    {

        RmgTimer *RT1 = new RmgTimer("VhMg");
        get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, ct.boundaryflag);

        delete(RT1);
    }


}
