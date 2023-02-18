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
        if( ct.force_grad_order != 0)  
        {
            ct.force_grad_order = 12;
        }
    }

}
