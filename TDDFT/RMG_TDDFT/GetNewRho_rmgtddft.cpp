/************************** SVN Revision Information **************************
 **    $Id: get_new_rho_local.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "RmgParallelFft.h"
#include "RmgGemm.h"
#include "blas_driver.h"



void GetNewRho_rmgtddft (double *psi, double *xpsi, double *rho, double *rho_matrix, int numst)
{
    int idx;

    /* for parallel libraries */

    int st1;
    double *rho_temp;

    double one = 1.0, zero = 0.0;
    int pbasis = get_P0_BASIS();

    rho_temp = new double[pbasis];

    for(idx = 0; idx < pbasis; idx++)rho_temp[idx] = 0.0;

    if(numst > 0)
    {
        

        RmgGemm ("N", "N", pbasis, numst, numst, one, 
                psi, pbasis, rho_matrix, numst, zero, xpsi, pbasis);

        my_sync_device();
        for(st1 = 0; st1 < numst; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += psi[st1 * pbasis + idx] * xpsi[st1 * pbasis + idx];

    }


    /* Interpolate onto fine grid, result will be stored in rho*/
    switch (ct.interp_flag)
    {
        case CUBIC_POLYNOMIAL_INTERPOLATION:
            pack_rho_ctof (rho_temp, rho);
            break;
        case PROLONG_INTERPOLATION:
            mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;
        case FFT_INTERPOLATION:
            FftInterpolation (*Rmg_G, rho_temp, rho, Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
     //       printf("\n Fftint not yet \n");
     //       exit(0);
            break;

        default:

            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }


    if(!ct.norm_conserving_pp) {

        printf("\n tddft not programed for ultrasoft \n");
        exit(0);
    }

    /* Check total charge. */
    ct.tcharge = ZERO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.img_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrpancy because of interpolation */
    double t1 = ct.nel / ct.tcharge;
    rmg_printf ("normalization constant-1 for new charge is %f\n", t1-1);
    for(int i = 0;i < FP0_BASIS;i++) rho[i] *= t1;

    delete [] rho_temp;


}
