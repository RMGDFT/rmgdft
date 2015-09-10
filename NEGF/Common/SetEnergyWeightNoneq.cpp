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


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"




/* This function returns a pointer to a block of memory of size nelem. */
void SetEnergyWeightNoneq (std::complex<double> *ene, std::complex<double> *weight, double EF1,
                           double EF2, int *nenergy)
{

    double a, b, tem;
    std::complex<double> distri1, distri2, ctem, II, ene_tem;
    int i, nen;
    int nmax_gq2;
    double KT, DELTA;
    double *xc, *wc;

    II.real(0.0);
    II.imag(1.0);
    nmax_gq2 = cei.nmax_gq2;
    KT = cei.KT;
    DELTA = cei.DELTA;

    xc = new double[nmax_gq2];
    wc = new double[nmax_gq2];

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

        ene_tem = xc[i] + II * DELTA;


        ctem = (ene_tem - EF1) / KT;
        distri1 = 1.0/( 1.0 + exp (ctem));

        ctem = (ene_tem - EF2) / KT;
        distri2 = 1.0/( 1.0 + exp (ctem));

        tem = abs(distri1 - distri2) ;
        if (tem > 0.000001)
        {

            ene[nen] = xc[i] + II * DELTA;

            weight[nen] = (distri2 - distri1) * wc[i]/2.0;
            nen++;
        }
    }

    *nenergy = nen;

    delete [] xc;
    delete [] wc;

    if (pct.gridpe == 0)
    {
        printf ("\n set_energy_weigh_ne done %d", *nenergy);
        printf ("\n eneR   eneI   weightR   weightI ");
        for (nen = 0; nen < *nenergy; nen++)
            printf ("\n  %f %f %f %f ", real(ene[nen]), imag(ene[nen]),
                    real(weight[nen]), imag(weight[nen]));
    }

}
