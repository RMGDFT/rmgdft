/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


void lcao_get_rho (double * arho_f)
{

    int ion, idx, FP0_BASIS;
    double t1, t2, difference;

    FP0_BASIS = get_FP0_BASIS();
    
    /* Initialize the compensating charge array and the core charge array */
    for (idx = 0; idx < get_FP0_BASIS(); idx++) arho_f[idx] = 0.0;

    // Accumulate atomic charge
    for (ion = 0; ion < pct.num_loc_ions; ion++) {
        for(idx = 0;idx < FP0_BASIS;idx++) arho_f[idx] += pct.localatomicrho[ion * FP0_BASIS + idx];
    }

    /* Check total charge. */
    t2 = 0.0;

    for (idx = 0; idx < get_FP0_BASIS(); idx++) t2 += arho_f[idx];

    t2 = get_vel_f() *  real_sum_all (t2, pct.grid_comm);

    t1 = ct.nel / t2;
    
    difference = fabs(t1 - 1.0);
    
    if ((ct.verbose == 1) || (difference > 0.05))
    {
	if (pct.imgpe == 0)
	    printf ("\n LCAO initializationNormalization constant for initial atomic charge is %f\n", t1);
    }

    for(idx = 0;idx < FP0_BASIS;idx++) arho_f[idx] *= t1;

}

