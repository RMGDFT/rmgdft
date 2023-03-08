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
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "GlobalSums.h"
#include "RmgException.h"
#include "InputKey.h"
#include "InputOpts.h"

void CheckSetDefault(void)
{
    if(ct.cell_relax)
    {
        ct.stress = true;

        bool cell_move_setup = false;
        for (int i = 0; i < 9; i++)
        {
            if(ct.cell_movable[i] != 0) cell_move_setup = true;
        } 

        if(!cell_move_setup)
        {
            for (int i = 0; i < 3; i++)
            {
                if(abs(Rmg_L.a0[i]) > 1.0e-5) ct.cell_movable[0*3 + i] = 1;
                if(abs(Rmg_L.a1[i]) > 1.0e-5) ct.cell_movable[1*3 + i] = 1;
                if(abs(Rmg_L.a2[i]) > 1.0e-5) ct.cell_movable[2*3 + i] = 1;
            }
        }

    }
    
    if(ct.stress)
    {
        ct.kohn_sham_fd_order = 12;
        ct.force_grad_order = 0;
        if( ct.force_grad_order != 0)  
        {
            ct.force_grad_order = 12;
        }
    }

    if(ct.spinorbit)
    {
        bool pp_has_so = false;
        for(int isp = 0;isp < ct.num_species;isp++)
        {
            SPECIES *sp = &Species[isp];
            if(sp->is_spinorb) pp_has_so = true;
        }
        if(!pp_has_so)
        {
            rmg_error_handler (__FILE__, __LINE__, "no pseudopotential has spin-orbit.\n");
        }
        ct.noncoll = true;
    }
    if(ct.noncoll) ct.is_ddd_non_diagonal = true;

    ct.norm_conserving_pp = true;
    int nc_count = 0;
    int us_count = 0;
    for(int isp = 0;isp < ct.num_species;isp++)
    {
        SPECIES *sp = &Species[isp];
        if(sp->is_norm_conserving)
        {
            nc_count++;
        }
        else
        {
            us_count++;
            ct.norm_conserving_pp = false;
        }
    }
    if(nc_count && us_count)
    {
        rmg_error_handler (__FILE__, __LINE__, "Mixing norm conserving and ultrasoft pseudopotentials is not supported. Check your input files.\n");
    }

    // For USPP force a minimum of 2
    if(!ct.norm_conserving_pp) ct.FG_RATIO = std::max(2, ct.FG_RATIO);

       /* Initialize symmetry stuff */

    ct.is_gamma = true;
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[0] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[1] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[2] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[0] == 0);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[1] == 0);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[2] == 0);
    ct.is_gamma = ct.is_gamma && (!ct.noncoll);
    if(ct.is_use_symmetry == 2)
    {
        ct.is_use_symmetry = 1;
        if(ct.is_gamma) ct.is_use_symmetry = 0;
    }
    if(!ct.is_gamma && ct.is_use_symmetry == 0 && pct.worldrank == 0)
    {
        printf("\n **********************************************************************");
        printf("\n WARNING: You turned off the symmetry for non-gamma point calculation ");
        printf("\n          use_symmetry = \"true\" to turn on symmetry");
        printf("\n *********************************************************************\n\n");
    }


}
