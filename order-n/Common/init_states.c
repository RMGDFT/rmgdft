/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/init_states.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2005  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 *	STATE *init_states()
 * 		setup state information
 *    
 * INPUTS
 *   nothing
 * OUTPUT
 *   pointer to STATE structure
 * PARENTS
 *   main.c
 * CHILDREN
 * 
 * SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "main.h"

#define MAX_NOCC 10


STATE *init_states ()
{
    int ii, is, ns, ik, nk;
    char *tbuf;
    int count_states = 0;

    struct
    {
        int n;
        double occ;
    } occ[MAX_NOCC];

    int nocc = 0;
    int repeat_occ = (strcmp (ct.occupation_str, "") != 0);


    /* calculate total number of electrons in the system from pseudopotential information */
    ct.ionic_charge = 0.0;

    for (ii = 0; ii < ct.num_ions; ii++)
        ct.ionic_charge += ct.sp[ct.ions[ii].species].zvalence;


    nk = ct.num_kpts;

    /* get repeat count occupation info */
    if (repeat_occ)
    {
        int n;

        ct.nel = 0;
        tbuf = ct.occupation_str;

        while ((n = strtol (tbuf, &tbuf, 10)) > 0)
        {
            count_states += n;

            if (nocc == MAX_NOCC)
                error_handler ("Too many blocks in repeat count for state occupations");

            occ[nocc].n = n;
            occ[nocc].occ = strtod (tbuf, &tbuf);

            ct.nel += n * occ[nocc].occ;

            nocc++;
        }


        /* ct.num_states is zero, if it has not been read from control file */
        if (ct.num_states)
        {
            if (ct.num_states != count_states)
            {
                printf ("\n ct.num_states:%d count_states:%d", ct.num_states, count_states);
                error_handler
                    ("states_count_and_occupation does not match specified number of states");
            }
        }
        else
            ct.num_states = count_states;


        /* calculate the compensating background charge
           for charged supercell calculations */

        ct.background_charge = ct.nel - ct.ionic_charge;

    }                           /*end if (repeat_occ) */
    else
    {
        ct.nel = ct.ionic_charge + ct.background_charge;
    }


    /* set the initial occupations of the states */

    if (repeat_occ)
    {
        int i;
        ns = 0;

        for (i = 0; i < nocc; i++)
        {
            for (is = ns; is < ns + occ[i].n; is++)
            {
                for (ik = 0; ik < nk; ik++)
                    states[ik * ct.num_states + is].occupation = occ[i].occ;

            }
            ns += occ[i].n;
        }
    }

    else
    {
        double oc;
        double ne = ct.nel;

        for (is = 0; is < ct.num_states; is++)
        {

            oc = 0.0;

            if (ne >= 2.0)
                oc = 2.0;
            else if (ne >= 0.0)
                oc = ne;

            ne -= oc;

            for (ik = 0; ik < nk; ik++)
                states[ik * ct.num_states + is].occupation = oc;

        }
    }


    if (pct.thispe == 0)
    {
        printf ("\n");
        printf ("total pseudopotential charge =  %8.3f e\n", ct.ionic_charge);
        printf ("total electronic charge      =  %8.3f e\n", -ct.nel);
        printf ("total system charge          =  %8.3f e\n", -ct.background_charge);
    }


}
