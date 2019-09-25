/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>



#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"

#include "transition.h"
#include "prototypes_on.h"

#include "blas.h"
#include "RmgParallelFft.h"
#include "BaseGrid.h"

void GetNewRho_proj(LocalObject<double> &Phi, LocalObject<double> &HPhi, double *rho, double *rho_matrix_local)
{

    double *rho_temp;

    double one = 1.0, zero = 0.0;

    RmgTimer *RT0 = new RmgTimer("3-get_new_rho");

    int pbasis = Rmg_G->get_P0_BASIS(1);
    rho_temp = new double[pbasis];

        
    for(int idx = 0; idx < pbasis; idx++)rho_temp[idx] = 0.0;
    if(Phi.num_thispe > 0)
    {

        int num_orb = Phi.num_thispe;
        dgemm ("N", "N", &pbasis, &num_orb, &num_orb, &one, 
                Phi.storage_proj, &pbasis, rho_matrix_local, &num_orb,
                &zero, HPhi.storage_proj, &pbasis);


        for(int st1 = 0; st1 < num_orb; st1++)
            for(int idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += Phi.storage_proj[st1*pbasis + idx] * HPhi.storage_proj[st1 * pbasis + idx];
    }



    RmgTimer *RT4 = new RmgTimer("3-get_new_rho: interpolation");
    /* Interpolate onto fine grid, result will be stored in rho*/
    switch (ct.interp_flag)
    {
        case PROLONG_INTERPOLATION:
            mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;
        case FFT_INTERPOLATION:
            //FftFilter(work, *coarse_pwaves, ct.cparm, LOW_PASS);  // limit to G-vectors within the inscribed sphere
            FftInterpolation (*Rmg_G, rho_temp, rho, Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
            break;
        default:
            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            printf("\n ct.interp_flag = %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }

    delete [] rho_temp;

    delete(RT4);

    RmgTimer *RT5 = new RmgTimer("3-get_new_rho: augmented");

//    RhoAugmented(rho, rho_matrix_local);


    delete(RT5);

    delete(RT0);


}
