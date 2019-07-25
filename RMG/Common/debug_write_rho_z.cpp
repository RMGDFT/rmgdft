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


#include <stdlib.h>
#include <stdio.h>
#include "common_prototypes.h"
#include "main.h"

void debug_write_rho_z (double * rhoz)
{
    int k;
    int basis;
    FILE *ftpr;

    if (pct.gridpe == 0)
    {
        my_fopen (ftpr, "rho_soft_separate.txt", "a+");

        basis = (get_PX0_GRID() / 2) * get_PY0_GRID() * get_PZ0_GRID() + (get_PY0_GRID() / 2) * get_PZ0_GRID();
        for (k = 0; k < get_PZ0_GRID(); k++)
            fprintf (ftpr, "%d   %f \n", k, rhoz[basis + k]);

        fclose (ftpr);
    }
}
