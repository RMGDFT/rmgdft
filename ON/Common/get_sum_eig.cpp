/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                        get_sum_eig.c



*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


double get_te_ion_ion()
{
    int i, j;
    double t1, energy = 0., r;
    ION *iptr1, *iptr2;



    /* Evaluate total ion-ion energy */
    energy = 0.;
    for (i = 0; i < ct.num_ions; i++)
    {

        energy -= (Species[Atoms[i].species].zvalence *
                   Species[Atoms[i].species].zvalence / Species[Atoms[i].species].rc);

    }                           /* end for */
    energy /= sqrt(2. * PI);


    for (i = 0; i < ct.num_ions; i++)
    {

        iptr1 = &Atoms[i];
        for (j = i + 1; j < ct.num_ions; j++)
        {

            iptr2 = &Atoms[j];

            r = minimage1(iptr1->crds, iptr2->crds);

            t1 = sqrt(Species[iptr1->species].rc * Species[iptr1->species].rc +
                      Species[iptr2->species].rc * Species[iptr2->species].rc);

            energy += Species[iptr1->species].zvalence *
                Species[iptr2->species].zvalence * erfc(r / t1) / r;


        }                       /* end for */

    }                           /* end for */

    return energy;
}

/***********************************************************/

double get_sum_eig(STATE * states)
{
    double eigsum = 0.0, t1;
    int kpt, state;

    /* Loop over states and get sum of the eigenvalues */

    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        t1 = ZERO;
        for (state = 0; state < ct.num_states; state++)
        {

            t1 += states[state + kpt * ct.num_states].occupation[0] *
                states[state + kpt * ct.num_states].eig[0];

        }                       /* end for */
        eigsum += t1 * ct.kp[kpt].kweight;
    }

    return eigsum;
}


/***********************************************************/

double get_Exc(double *rho, double *rhocore)
{
    double *exc, *nrho, esum;
    int idx;


    /* Grab some memory */
    /* begin shuchun wang */
    my_malloc_init( exc, ct.vh_pbasis, double );
    my_malloc_init( nrho, ct.vh_pbasis, double );

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        nrho[idx] = rhocore[idx] + rho[idx];

    get_vxc_exc (nrho,  nrho, vxc, exc, ct.xctype);


    esum = 0.;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {

        esum += (rhocore[idx] + rho[idx]) * exc[idx];

    }                           /* end for */

    /* Release our memory */
    my_free(nrho);
    my_free(exc);


    esum = get_vel_f() * real_sum_all(esum, pct.grid_comm);
    /* end shuchun wang */

    return esum;
}


/***********************************************************/
