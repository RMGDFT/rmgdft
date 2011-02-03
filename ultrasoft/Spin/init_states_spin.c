/************************** SVN Revision Information **************************
 **    $Id: init_states.c 1066 2009-08-31 18:41:09Z froze $    **
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


STATE *init_states_spin ()
{
    int ii, is, ns, ik, nk, temp;
    STATE *states;
    char *tbuf;
    int count_states_up = 0, count_states = 0;
    int count_states_down=0;

    struct
    {
        int n;
        double occ;
    } occ_up[MAX_NOCC], occ_down[MAX_NOCC];


    int nocc_up = 0, nocc_down = 0;
    int repeat_occ ;    

    /* calculate total number of electrons in the system from pseudopotential information */
    ct.ionic_charge = 0.0;

    for (ii = 0; ii < ct.num_ions; ii++)
        ct.ionic_charge += ct.sp[ct.ions[ii].species].zvalence;


    repeat_occ =( (strcmp(ct.occupation_str_spin_up, "") != 0) && (strcmp(ct.occupation_str_spin_down, "")!= 0) );    

    /* get repeat count occupation info */
    if (repeat_occ)
    {
        int n;

        ct.nel = 0;
        
	/* count the fixed occupations for spin up states*/
        tbuf = ct.occupation_str_spin_up;
        while ((n = strtol (tbuf, &tbuf, 10)) > 0)
        {
            count_states_up += n;

            if (nocc_up == MAX_NOCC)
                error_handler ("Too many blocks in repeat count for state occupations");  /* two block example  3 1.0  4 0.0*/

            occ_up[nocc_up].n = n;
            occ_up[nocc_up].occ = strtod (tbuf, &tbuf);

            ct.nel += n * occ_up[nocc_up].occ;

            nocc_up++;
        }

        /* count the fixed occupations for spin down states*/ 
        tbuf = ct.occupation_str_spin_down;
        while ((n = strtol (tbuf, &tbuf, 10)) > 0)
        {
            count_states_down += n;

            if (nocc_down == MAX_NOCC)
                error_handler ("Too many blocks in repeat count for state occupations");

            occ_down[nocc_down].n = n;
            occ_down[nocc_down].occ = strtod (tbuf, &tbuf);

            ct.nel += n * occ_down[nocc_down].occ;

            nocc_down++;
        }

	/* ct.num_states_up is not zero, if it has been read from control file */
        if (ct.num_states_up)
        {	
            if ( ct.num_states_up != count_states_up)
            {
            	printf ("\n ct.num_states_up:%d count_states_up:%d", ct.num_states_up, count_states_down);
                error_handler("states_count_and_occupation does not match specified number of states");
            }
        }
 	else 
	    ct.num_states_up=count_states_up;       
                
	/* ct.num_states_down is zero, if it has been read from control file */
        if (ct.num_states_down)
        {	
            if ( ct.num_states_up != count_states_up)
            {
                printf ("\n ct.num_states_up:%d count_states_up:%d", ct.num_states_up, count_states_down);
                error_handler("states_count_and_occupation does not match specified number of states");
            }
        }
	else
	    ct.num_states_down=count_states_down;


        /* calculate the compensating background charge
           for charged supercell calculations */

        ct.background_charge = ct.nel - ct.ionic_charge;

    }                           /*end if (repeat_occ) */
    else
    {
        ct.nel = ct.ionic_charge + ct.background_charge;
        ct.num_states_up = (int) ceil(ct.nel / 2) + ct.num_unocc_states;
        ct.num_states_down = (int)  ceil(ct.nel / 2) + ct.num_unocc_states;
    }


    /* Allocate memory for the states */
    if (verify ("calculation_mode", "Band Structure Only"))
        nk = 1;
    else
        nk = ct.num_kpts;


    temp = (ct.num_states_up >= ct.num_states_down) ? ct.num_states_up : ct.num_states_down;

    my_malloc (states, temp * nk, STATE);


    /* set the initial occupations of the states */

    if (repeat_occ)
    {
        int i;
        
	/* the rank 0 in spin communicator will handle spin up */    
        if (pct.thisspin==0)
	{
		ns = 0;
        	for (i = 0; i < nocc_up; i++)
        	{
            		for (is = ns; is < ns + occ_up[i].n; is++)
            		{
                		for (ik = 0; ik < nk; ik++)
                    		states[ik * ct.num_states_up + is].occupation = occ_up[i].occ;

            		}
            		ns += occ_up[i].n;
        	}

		ns = 0;
        	for (i = 0; i < nocc_down; i++)
        	{
            		for (is = ns; is < ns + occ_down[i].n; is++)
            		{
                		for (ik = 0; ik < nk; ik++)
                    		states[ik * ct.num_states_down + is].occupation_oppo = occ_down[i].occ;

            		}
            		ns += occ_down[i].n;
        	}
		ct.num_states=ct.num_states_up;
		ct.num_states_oppo=ct.num_states_down;
    	}

	/* the rank 1 in spin communicator will handle spin down */    
        else if (pct.thisspin==1)
	{
		ns = 0;
        	for (i = 0; i < nocc_down; i++)
        	{
            		for (is = ns; is < ns + occ_down[i].n; is++)
            		{
                		for (ik = 0; ik < nk; ik++)
                    		states[ik * ct.num_states_down + is].occupation = occ_down[i].occ;
            		}
            		ns += occ_down[i].n;
        	}

		ns = 0;
        	for (i = 0; i < nocc_up; i++)
        	{
            		for (is = ns; is < ns + occ_up[i].n; is++)
            		{
                		for (ik = 0; ik < nk; ik++)
                    		states[ik * ct.num_states_up + is].occupation_oppo = occ_up[i].occ;
            		}
            		ns += occ_up[i].n;
        	}
		ct.num_states=ct.num_states_down;
		ct.num_states_oppo=ct.num_states_up;
    	}

    }
    
    else   /* occupation is not specified by input file*/
    {
        double oc;
	double ne_up, ne_down, ne, ne_oppo ;

	ne_up = ct.nel/2.0;
	ne_down = ct.nel - ne_up;
       	
        
	/* set rank 0 in spin communicator to handle spin up*/
	if (pct.thisspin ==0)
	{
		ct.num_states=ct.num_states_up;
		ct.num_states_oppo=ct.num_states_down;
		ne = ne_up;
		ne_oppo = ne_down;
	}
	/* set rank 0 in spin communicator to handle spin up*/
	else if (pct.thisspin ==1)
	{
		ct.num_states=ct.num_states_down;
		ct.num_states_oppo=ct.num_states_up;
		ne = ne_down;
		ne_oppo = ne_up;
	}

        /* initialize occupations for the spin the processor handles*/
        for (is = 0; is < ct.num_states; is++)
        {
            oc = 0.0;
            if (ne >= 1.0)
                oc = 1.0;
            else if (ne >= 0.0)
                oc = ne;
            ne -= oc;

            for (ik = 0; ik < nk; ik++)
                states[ik * ct.num_states + is].occupation = oc;
        }

        /* initialize occupations for the opposite spin*/
        for (is = 0; is < ct.num_states_oppo; is++)
        {
            oc = 0.0;
            if (ne_oppo >= 1.0)
                oc = 1.0;
            else if (ne_oppo >= 0.0)
                oc = ne_oppo;
            ne_oppo -= oc;

            for (ik = 0; ik < nk; ik++)
                states[ik * ct.num_states + is].occupation_oppo = oc;
        }
    }

    if (pct.imgpe == 0 )
    {
        printf ("\n");
        printf ("total pseudopotential charge =  %8.3f e\n", ct.ionic_charge);
        printf ("total electronic charge      =  %8.3f e\n", -ct.nel);
        printf ("total system charge          =  %8.3f e\n", -ct.background_charge);
    }

    return states;
}

