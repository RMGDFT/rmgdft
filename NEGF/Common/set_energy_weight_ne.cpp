#include "negf_prototypes.h"
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


#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* This function returns a pointer to a block of memory of size nelem. */
void set_energy_weight_ne (std::complex<double> * ene, std::complex<double> * weight, double EF1,
                           double EF2, int *nenergy)
{

    double a, b;
    std::complex<double>  distri1, distri2, ctem;
    std::complex<double> I(0.0, 1.0);
    int i, nen;
    int nmax_gq2;
    double KT, DELTA;
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
        rmg_printf ("\n fff %f %f %d", a, b, nmax_gq2);

/* then determine the Gauss-Legand parameters */
    gauleg (a, b, xc, wc, nmax_gq2);

/* calculate the energy points and corresponding weight  */

    nen = 0;
    for (i = 0; i < nmax_gq2; i++)
    {
        ctem = xc[i] + I * DELTA;
        distri_fermi (ctem, EF1, &distri1);
        distri_fermi (ctem, EF2, &distri2);
        {

            ene[nen] = xc[i] + I * DELTA;

            weight[nen] = (distri2 - distri1) * wc[i]/2.0;
            nen++;
        }
    }

    *nenergy = nen;

    my_free(xc);
    my_free(wc);
    if (pct.gridpe == 0)
    {
        rmg_printf ("\n set_energy_weigh_ne done %d", *nenergy);
        rmg_printf ("\n eneR   eneI   weightR   weightI ");
        for (nen = 0; nen < *nenergy; nen++)
            rmg_printf ("\n  %f %f %f %f ", std::real(ene[nen]), std::imag(ene[nen]),
                    std::real(weight[nen]), std::imag(weight[nen]));
    }

}
