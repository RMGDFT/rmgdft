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


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "transition.h"


// Returns the total number of atomic orbitals including the m-dependence.
// Used during LCAO intialization in order to determine required memory allocations.
int CountAtomicOrbitals(void)
{

    int total_atomic_orbitals = 0;
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        ION *iptr = &Atoms[ion];

        /* Get species type */
        SPECIES *sp = &Species[iptr->species];

        for (int ip = 0; ip < sp->num_atomic_waves; ip++) {
            if(sp->atomic_wave_oc[ip] > 0.0) {
                int l = sp->atomic_wave_l[ip];
                double jj = sp->atomic_wave_j[ip];

                if( sp->is_spinorb )
                   total_atomic_orbitals += (int)(2 * jj + 1 + 0.1);
                else if(ct.noncoll)
                    total_atomic_orbitals += 2*(2 * l + 1);
                else
                    total_atomic_orbitals += 2 * l + 1;
            }
        }

    }

    return total_atomic_orbitals;
}
