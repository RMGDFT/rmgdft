/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_milliken.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void get_milliken(STATE *states)
 *   Computes milliken population densities.
 * INPUTS
 *   states: point to orbital structure (see main.h)
 * OUTPUT
 *   print out milliken population
 * PARENTS
 *   fastrlx.c main.c
 * CHILDREN
 *   app_nl.c gather_psi.c get_nlop.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"




void get_milliken (STATE * states)
{

    int istate, pbasis, ion, ip;
    int channel;
    REAL *tmp_psi, *work;
    STATE *sp;
    ION *iptr;
    REAL charge[MAX_IONS];      /* mullikan charges */


    pbasis = (PX0_GRID * PY0_GRID * PZ0_GRID);
    my_malloc (tmp_psi, 2 * pbasis, REAL);
    my_malloc (work, 2 * pbasis, REAL);


    /* Zero the mullikan charges */
    for (ion = 0; ion < ct.num_ions; ion++)
        charge[ion] = 0.0;

    /* Generate the projectors */
    get_nlop (0);

    /* Now loop over states and apply the milliken projectors to all the
     * states. */

    for (istate = 0; istate < ct.num_states; istate++)
    {

        sp = &states[istate];
        gather_psi (tmp_psi, NULL, sp, 0);
        app_nl_eig (tmp_psi, NULL, work, NULL, sp->istate, FALSE, sp->kidx, 0);
        if (ct.milliken == 2)
        {
            if (pct.thispe == 0)
            {

                printf ("\n\nMilliken Population Information\n");

                printf ("Energy of state %d: %7.2f; Occupation: %4.2f\n",
                        istate, sp->eig * Ha_eV, sp->occupation);
            }                   /* if(pct.thispe == 0) */
        }                       /* if (ct.milliken==2 */
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Generate ion pointer */
            iptr = &ct.ions[ion];

            for (ip = 0; ip < pct.prj_per_ion[ion]; ip++)
            {

                if (ip == 0)
                {
                    channel = 0;
                }
                else if (ip < 4)
                {
                    channel = 1;
                }
                else
                {
                    channel = 2;
                }
                if (ct.milliken == 2)
                {
                    if (pct.thispe == 0)
                    {

                        printf ("State %d, ion %d, species %d, angular momentum %d, %12.6f\n",
                                istate, ion, iptr->species + 1, channel,
                                iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                               istate * ct.max_nl + ip]);
                    }           /* if(pct.thispe == 0) */
                }               /* if (ct.milliken==2) */

                charge[ion] +=
                    iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                   istate * ct.max_nl +
                                   ip] * iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states *
                                                        ct.max_nl + istate * ct.max_nl +
                                                        ip] * sp->occupation;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    if (ct.milliken > 0)
    {
        if (pct.thispe == 0)
        {
            printf ("\n\nMilliken Charges\n");
            for (ion = 0; ion < ct.num_ions; ion++)
            {
                /* Generate ion pointer */
                iptr = &ct.ions[ion];

                printf ("Ion %d, species %d, charge = %f\n", ion, iptr->species + 1, charge[ion]);
            }                   /* end for */
        }                       /* end if(pct.thispe == 0) */
    }                           /* if (ct.milliken>0) */
    ct.milliken = FALSE;

    my_free (work);
    my_free (tmp_psi);

}                               /* end get_milliken */

/******/
