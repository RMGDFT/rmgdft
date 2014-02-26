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

        energy -= (ct.sp[ct.ions[i].species].zvalence *
                   ct.sp[ct.ions[i].species].zvalence / ct.sp[ct.ions[i].species].rc);

    }                           /* end for */
    energy /= sqrt(2. * PI);


    for (i = 0; i < ct.num_ions; i++)
    {

        iptr1 = &ct.ions[i];
        for (j = i + 1; j < ct.num_ions; j++)
        {

            iptr2 = &ct.ions[j];

            r = minimage1(iptr1->crds, iptr2->crds);

            t1 = sqrt(ct.sp[iptr1->species].rc * ct.sp[iptr1->species].rc +
                      ct.sp[iptr2->species].rc * ct.sp[iptr2->species].rc);

            energy += ct.sp[iptr1->species].zvalence *
                ct.sp[iptr2->species].zvalence * erfc(r / t1) / r;


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
    double *work, *nrho, esum;
    int idx;

    rmg_double_t tem, tem1, tem2;

    /* Grab some memory */
    /* begin shuchun wang */
    my_malloc_init( work, ct.vh_pbasis, rmg_double_t );
    my_malloc_init( nrho, ct.vh_pbasis, rmg_double_t );

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        nrho[idx] = rhocore[idx] + rho[idx];


    /* Evaluate XC energy correction terms */
    switch (ct.xctype)
    {
    case LDA_PZ81:             /* LDA Perdew Zunger 81 */
        exclda_pz81(nrho, work, get_FP0_BASIS());
        break;

    case GGA_BLYP:             /* GGA X-Becke C-Lee Yang Parr */
        xcgga(nrho, vxc, work, ct.xctype);
        break;

    case GGA_XB_CP:            /* GGA X-Becke C-Perdew */
        xcgga(nrho, vxc, work, ct.xctype);
        break;

    case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */

        xcgga(nrho, vxc, work, ct.xctype);
        break;

    case GGA_PBE:              /* GGA Perdew, Burke, Ernzerhof */

        xcgga(nrho, vxc, work, ct.xctype);
        break;

    default:
        error_handler("Unknown exchange-correlation functional");

    }                           /* end switch */


    esum = 0.;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {

        esum += (rhocore[idx] + rho[idx]) * work[idx];

    }                           /* end for */

    /* Release our memory */
    my_free(nrho);
    my_free(work);


    esum = get_vel_f() * real_sum_all(esum, pct.grid_comm);
    /* end shuchun wang */

    return esum;
}


/***********************************************************/
