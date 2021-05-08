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
    int ii, is, ns, ik, nk, idx, i, j;
    STATE *states;
    int count_states[2]={0,0}, nocc[2]={0,0}, num_states_spf[2], nspin = (ct.spin_flag + 1);
    char *tbuf[2];

    struct
    {
        int n;
        double occ;
    } occ[nspin * MAX_NOCC];

    int repeat_occ;


    /* calculate total number of electrons in the system from pseudopotential information */
    ct.ionic_charge = 0.0;

    for (ii = 0; ii < ct.num_ions; ii++)
        ct.ionic_charge += Species[Atoms[ii].species].zvalence; 


    if (ct.spin_flag)
    {
    	repeat_occ =( (strcmp(ct.occupation_str_spin_up, "") != 0) && (strcmp(ct.occupation_str_spin_down, "")!= 0) );
    	if( (strcmp(ct.occupation_str_spin_up, "") != 0) + (strcmp(ct.occupation_str_spin_down, "")!= 0) == 1 )
		error_handler ("Fixed occupation for both spin up and down must be specified !!!");
	tbuf[0] = ct.occupation_str_spin_up;
	tbuf[1] = ct.occupation_str_spin_down;
	num_states_spf[0] = 0;
	num_states_spf[1] = 0;
    }	
    else
    { 
        repeat_occ = (strcmp (ct.occupation_str, "") != 0);
	num_states_spf[0] = 0;
	tbuf[0] = ct.occupation_str;
    }


    /* get repeat count occupation info from control file*/
    if (repeat_occ)
    {
        int n;
       	ct.nel = 0; 
		
        for (idx = 0; idx < nspin; idx++) 
	{
		/* count the fixed occupations of states */
        	while ((n = strtol (tbuf[idx], &tbuf[idx], 10)) > 0)
        	{
            		count_states[idx] += n;
            		if (nocc[idx] == MAX_NOCC)
                		error_handler ("Too many blocks in repeat count for state occupations");  
				/* two block example  3 1.0  4 0.0*/
            		occ[nocc[idx] + MAX_NOCC * idx].n = n;
            		occ[nocc[idx] + MAX_NOCC * idx].occ = strtod (tbuf[idx], &tbuf[idx]);
            		ct.nel += n * occ[nocc[idx] + MAX_NOCC * idx].occ;
            		nocc[idx]++;
        	}
	       
		num_states_spf[idx] = count_states[idx]; 
	} 
    	
        if ( (nspin == 2) && (num_states_spf[0] != num_states_spf[1]) )
	{       
		printf("number of states for spin up: %d, number of states for spin down %d\n", num_states_spf[0], num_states_spf[1]);
		error_handler("num_of_states_spin_up not equal to num_states_spin_down, you are wasting memory address for extra STATE structures !");
	}
        
	ct.background_charge = ct.nel - ct.ionic_charge; 

    } 
    else     /* in case no fixed occupations available, calculate number of states */
    {
        ct.nel = ct.ionic_charge + ct.background_charge;
        for (idx = 0; idx < nspin; idx++)
        	num_states_spf[idx] = (int) ceil(0.5 * ct.nel) + ct.num_unocc_states;
    }
    


    /* re-assign the number of states for global variables */
    if(ct.num_states <= 0 ) ct.num_states = num_states_spf[0];


    if (nspin == 2)
    {
	ct.num_states_up = num_states_spf[0];
	ct.num_states_down = num_states_spf[1];
    }

    /* Allocate memory for the states */
    nk = ct.num_kpts; 

    states = new STATE[ct.num_states * nk];

    /* set the initial occupations of the states */
    if (repeat_occ)
    {
	        
	if (nspin ==1)
		pct.spinpe = 0;	
		
	for (idx = 0; idx < nspin; idx++)
	{
		ns = 0;  
		j = (idx + pct.spinpe) % 2;
		for (i = 0; i < nocc[j]; i++)
		{
			for (is = ns; is < ns + occ[i + j * MAX_NOCC].n; is++)
				for (ik = 0; ik < nk; ik++)
                    			states[ik * num_states_spf[j] + is].occupation[idx] = occ[i + j * MAX_NOCC].occ;
				
			ns += occ[i + j * MAX_NOCC].n;
		}
	}
    }

    else
    {       
	double ne[nspin], oc; 
	
	for (idx = 0; idx < nspin; idx++)
		ne[idx] = ct.nel / ((double) nspin);

	for (idx = 0; idx < nspin; idx++)
	{
		for (is = 0; is < ct.num_states; is++)
		{
			oc = 0.0;
			if ( ne[idx] >= (3.0 - nspin) )
				oc = (3.0 - nspin);
			else if (ne[idx] >= 0.0)
				oc = ne[idx];
			ne[idx] -= oc;
			for (ik = 0; ik < nk; ik++)
				states[ik * ct.num_states + is].occupation[idx] = oc;
		}
	}

    }


    /* Print out results to output file */ 
    printf ("\n");
    printf ("total pseudopotential charge =  %8.3f e\n", ct.ionic_charge);
    printf ("total electronic charge      =  %8.3f e\n", -ct.nel);
    printf ("total system charge          =  %8.3f e\n", -ct.background_charge);

    return states;
}
