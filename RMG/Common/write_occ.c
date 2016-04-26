/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "portability.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"


/* Writes occupations */
void write_occ (STATE * states)
{
    int i, idx, nspin = (ct.spin_flag + 1);

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

    	printf ("\n\n");

    }

}                               /* end write_occ */

/******/
