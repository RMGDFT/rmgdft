/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"



void force(REAL * rho, REAL * rhoc, REAL * vh, REAL * vxc, REAL * vnuc, STATE * states)
{
    int ion, st, kpt, idx;
    STATE *sp;

    REAL time1, time2;
    time1 = my_crtc();

    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    /* Zero out forces */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        ct.ions[ion].force[ct.fpt[0]][0] = ZERO;
        ct.ions[ion].force[ct.fpt[0]][1] = ZERO;
        ct.ions[ion].force[ct.fpt[0]][2] = ZERO;

    }                           /* end for */

    /* Get the ion-ion component and store. */
    iiforce();

    if (pct.gridpe == 0)
    {
        printf("\n iiforce");
        for (ion = 0; ion < ct.num_ions; ion++)
            printf("\n     %d  %10.7f %10.7f %10.7f",
                   ion, ct.ions[ion].force[0][0], ct.ions[ion].force[0][1],
                   ct.ions[ion].force[0][2]);
    }

    /* Add in the local */
    lforce(rho, vh);

    if (pct.gridpe == 0)
    {
        printf("\n lforce");
        for (ion = 0; ion < ct.num_ions; ion++)
            printf("\n     %d  %10.7f %10.7f %10.7f",
                   ion, ct.ions[ion].force[0][0], ct.ions[ion].force[0][1],
                   ct.ions[ion].force[0][2]);
    }


    /* Add in the non-local stuff */
    nlforce(vtot, states);

    if (pct.gridpe == 0)
    {
        printf("\n nlforce");
        for (ion = 0; ion < ct.num_ions; ion++)
            printf("\n     %d  %10.7f %10.7f %10.7f\n",
                   ion, ct.ions[ion].force[0][0], ct.ions[ion].force[0][1],
                   ct.ions[ion].force[0][2]);
    }



    /* The non-linear core correction part if any */
    nlccforce(rho, vxc);

#if !GAMMA_PT
    /* Now symmetrize the forces */
    if (!(ct.kp[0].kpt[0] == 0.0 && ct.kp[0].kpt[1] == 0.0 && ct.kp[0].kpt[2] == 0.0))
        symforce();
#endif

    time2 = my_crtc();
    rmg_timings(FORCE_TIME, time2 - time1, 0);

}


/******/
