/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/set_energy_weight.c ************
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
#include <complex>
#include <cmath>

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
void SetEnergyWeight (std::complex<double> * ene, std::complex<double> * weight, double EF, int *nenergy)
{

    int nen, nloop, i;
    double tem, emax, emin, rad, theta0, thetamax, thetamin;
    double center;
    std::complex<double>  ctem, distri, II;

    double KT, DELTA, DELTA2, GAMMA, EB;
    int size,ncircle, nmax_gq1;
    double *xc, *wc;

    II.real(0.0);
    II.imag(1.0);
    KT = cei.KT;
    EB = cei.EB;
    DELTA = cei.DELTA;
    GAMMA = cei.GAMMA;
    DELTA2 = cei.DELTA2;
    ncircle = cei.ncircle;
    nmax_gq1 = cei.nmax_gq1;

    size = std::max(ncircle, nmax_gq1); 
    xc = new double[size];
    wc = new double[size];


/*   first count how many poles in EF + i (2mu+1)PI *kt 
 */


    nen = 0;
    nloop = 10 * (int) (DELTA2 / (PI * KT));
    for (i = 0; i < nloop; i++)
    {
        tem = (2 * i + 1) * PI * KT;
        if (tem < DELTA2)
        {
            ene[nen] = EF + II * tem;
            weight[nen] = II * 2.0 * PI * KT;
            nen++;
        }
    }


/*  second, the (almost) semi-circle */

    emin = EB - EF;
    emax = -GAMMA + EF;
/*  radius of the circle  */
    rad = ((emax - emin) * (emax - emin) + DELTA2 * DELTA2) / (2.0 * (emax - emin));

    /*  opening angle of the semi-circle PI -- theta0  */
    center = emin + rad;        /* center of the semi-circle */
    theta0 = atan (DELTA2 / (emax - center));
    thetamax = PI;
    thetamin = theta0;
    gauleg (thetamin, thetamax, xc, wc, ncircle);

    for (i = 0; i < ncircle; i++)
    {
        theta0 = xc[i];
        ene[nen] = rad * exp ( II*theta0 ) + center + II * DELTA;

        ctem = (ene[nen] - EF) / KT;
        distri = 1.0/( 1.0 + exp (ctem));

        weight[nen] = II * distri * exp( II*theta0 ) * rad * wc[i];
        nen++;

    }

    /*  finally work on L part */
    emin = -GAMMA + EF;
    /*	emax = EF + GAMMA; */
    emax = EF + 2.5;
    gauleg (emin, emax, xc, wc, nmax_gq1);

    for (i = 0; i < nmax_gq1; i++)
    {

        ene[nen] = xc[i] + II * (DELTA2 + DELTA);

        ctem = (ene[nen] - EF) / KT;
        distri = 1.0/( 1.0 + exp (ctem));

        weight[nen] = -distri * wc[i];
        nen++;

    }

    *nenergy = nen;

    delete [] xc;
    delete [] wc;

    if (pct.gridpe == 0)
    {
        printf ("\n set_energy_weight done %d", *nenergy);
        printf ("\n    eneR   eneI   weightR   weightI");
        for (i = 0; i < nen; i++)
            printf ("\n     %f  %f %f  %f", real(ene[i]), imag(ene[i]), 
                    real(weight[i]), imag(weight[i]));
    }

}

/*******/
