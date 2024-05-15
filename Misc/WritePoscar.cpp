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


#include <stdio.h>
#include <time.h>
#include <math.h>

#include "common_prototypes.h"
#include "Functional.h"
#include "RmgParallelFft.h"
#include "transition.h"
#include "Functional.h"

/* Writes out header information */
void WritePoscar(FILE *FH, int md_step)
{
    if(md_step == 0)
    {
        fprintf(FH, "%s\n", ct.description.c_str());
        fprintf(FH, "1.0\n");
        fprintf(FH, "%f %f %f\n", Rmg_L.a0[0] * a0_A, Rmg_L.a0[1] * a0_A, Rmg_L.a0[2] * a0_A); 
        fprintf(FH, "%f %f %f\n", Rmg_L.a1[0] * a0_A, Rmg_L.a1[1] * a0_A, Rmg_L.a1[2] * a0_A); 
        fprintf(FH, "%f %f %f\n", Rmg_L.a2[0] * a0_A, Rmg_L.a2[1] * a0_A, Rmg_L.a2[2] * a0_A); 

        for (auto& sp : Species)
        {
            fprintf(FH, "%6.6s", sp.atomic_symbol);
        }
        fprintf(FH, "\n");

        for (auto& sp : Species)
        {
            fprintf(FH, "%6d  ", sp.num_atoms);
        }
        fprintf(FH, "\n");
    }
    fprintf(FH, "Cartesian configuration %d\n", md_step);

    for (int isp = 0; isp < (int)Species.size(); isp++)
    {
        for(auto& Atom : Atoms)
        {
            if(Atom.species == isp) 
            {
                fprintf(FH, "%16.8e %16.8e %16.8e\n", Atom.crds[0] * a0_A, Atom.crds[1] * a0_A, Atom.crds[2] * a0_A); 
            }
        }
    }

    fflush(FH);
}
