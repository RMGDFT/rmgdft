/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/set_energy_weight_ne.c ************
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 * FUNCTION
 *   void set_energy_weight(real eneR, real eneI, real weight, nenergy)
 *   set up the energies and weights and # of Green functions 
 * INPUTS
 *   
 * OUTPUT
 *   nothing
 * PARENTS
 *   too many
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* This function returns a pointer to a block of memory of size nelem. */
void set_energy_weight_ne (REAL * eneR, REAL * eneI, REAL * weightR, REAL * weightI, REAL EF1,
                           REAL EF2, int *nenergy)
{

    REAL a, b, distriR1, distriI1, distriR2, distriI2, tem;
    int i, nen;
    int nmax_gq2;
    REAL KT, DELTA;
    double *xc, *wc;

    nmax_gq2 = cei.nmax_gq2;
    KT = cei.KT;
    DELTA = cei.DELTA;

     my_malloc( xc, nmax_gq2, double );
     my_malloc( wc, nmax_gq2, double );

/*  First determine the energy region for the non-equilbrium part  */

    if (EF1 > EF2)
    {
        a = EF2 - 25.0 * KT;
        b = EF1 + 25.0 * KT;
    }
    else
    {
        a = EF1 - 25.0 * KT;
        b = EF2 + 25.0 * KT;
    }
    if (pct.gridpe == 0)
        printf ("\n fff %f %f %d", a, b, nmax_gq2);

/* then determine the Gauss-Legand parameters */
    gauleg (a, b, xc, wc, nmax_gq2);

/* calculate the energy points and corresponding weight  */

    nen = 0;
    for (i = 0; i < nmax_gq2; i++)
    {
        distri_fermi (xc[i], DELTA, EF1, &distriR1, &distriI1);
        distri_fermi (xc[i], DELTA, EF2, &distriR2, &distriI2);
        tem = (distriR1 - distriR2) * (distriR1 - distriR2);
        tem += (distriI1 - distriI2) * (distriI1 - distriI2);
        tem = sqrt (tem);
        if (tem > 0.000001)
        {

            eneR[nen] = xc[i];
            eneI[nen] = DELTA;

            weightR[nen] = (distriR2 - distriR1) * wc[i];
            weightI[nen] = (distriI2 - distriI1) * wc[i];
            nen++;
        }
    }

    *nenergy = nen;

    my_free(xc);
    my_free(wc);
    if (pct.gridpe == 0)
    {
        printf ("\n set_energy_weigh_ne done %d", *nenergy);
        printf ("\n eneR   eneI   weightR   weightI ");
        for (nen = 0; nen < *nenergy; nen++)
            printf ("\n  %f %f %f %f ", eneR[nen], eneI[nen], weightR[nen], weightI[nen]);
    }

}
