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

#include "main.h"



void write_pdb (void)
{
    int ion;
    ION *iptr;



    //printf("\n\n Updated PDB file");

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];


        printf ("\n");
        printf ("%-6s", iptr->pdb.record_name);
        printf ("%5d ", iptr->pdb.serial_num);
        printf ("%4s", iptr->pdb.name);
        printf ("%1s", iptr->pdb.altLoc);
        printf ("%3s ", iptr->pdb.resName);
        printf ("%1s", iptr->pdb.chainID);
        printf ("%4d", iptr->pdb.resSeq);
#if 1
        printf ("%1s   ", iptr->pdb.iCode);
        printf ("%8.3f", a0_A * iptr->crds[0]);
        printf ("%8.3f", a0_A * iptr->crds[1]);
        printf ("%8.3f", a0_A * iptr->crds[2]);
        //printf("%6s", iptr->pdb_occupancy);
        printf ("%6.2f", iptr->pdb.occupancy);
        //printf("%6s      ", iptr->pdb_tempFactor);
        printf ("%6.2f          ", iptr->pdb.tempFactor);
        //printf("%4s", iptr->pdb_segID);
        printf ("%2s", ct.sp[iptr->species].atomic_symbol);
        printf ("%2s", iptr->pdb.charge);
#endif

    }                           /*end for(ion = 0;ion < c->num_ions;ion++) */



    printf ("\nEND");
    printf ("\n\n");


}                               /* EOF */
