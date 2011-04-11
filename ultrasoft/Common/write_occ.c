/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_occ.c *****
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
 *   void write_occ(STATE *states)
 *   Writes occupations of eace state
 * INPUTS
 *   states: point to orbital structure
 * OUTPUT
 *   print out occupation numbers
 * PARENTS
 *   main.c quench.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"


/* Writes occupations */
void write_occ (STATE * states)
{
    int i, idx, nspin = (pct.spin_flag + 1);

    switch (ct.occ_flag)
    {
    case OCC_NONE:
        break;
    case OCC_FD:
        printf ("\nFERMI-DIRAC OCCUPATION WITH PARAMETERS:");
        printf ("\n  TEMP   = %14.8f", ct.occ_width);
        printf ("\n  MIXING = %14.8f", ct.occ_mix);
        break;
    case OCC_GS:
        printf ("\nGAUSSIAN OCCUPATION WITH PARAMETERS:");
        printf ("\n  TEMP   = %14.8f", ct.occ_width);
        printf ("\n  MIXING = %14.8f", ct.occ_mix);
        break;
    case OCC_EF:
        printf ("\nERROR_FUNCTION OCCUPATION WITH PARAMETERS:");
        printf ("\n  TEMP   = %14.8f", ct.occ_width);
        printf ("\n  MIXING = %14.8f", ct.occ_mix);
        break;
    default:
        error_handler ("unknown filling procedure");
    } 


    for (idx = 0; idx < nspin; idx++)
    {
	if (nspin == 1)
		printf ("\n\n  STATE OCCUPATIONS :\n");
	else if ((nspin == 2) && (idx == 0))
		printf ("\n\n  STATE OCCUPATIONS FOR SPIN UP:\n");
	else if ((nspin == 2) && (idx == 1))
		printf ("\n\n  STATE OCCUPATIONS FOR SPIN DOWN:\n"); 
    	
	for (i = 0; i < ct.num_states; i++)
        	printf (" %7.2f%s", states[i].occupation[idx], ((i % 10 == 9) ? "\n" : ""));

    	printf ("\n");

    }

}                               /* end write_occ */

/******/
