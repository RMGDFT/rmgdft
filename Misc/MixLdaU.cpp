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

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <unordered_map>
#include <typeinfo>
#include "const.h"
#include "InputKey.h"
#include "common_prototypes.h"
#include "RmgParallelFft.h"
#include "RmgException.h"
#include "transition.h"
#include "blas.h"
#include "GlobalSums.h"
#include "Functional.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <boost/circular_buffer.hpp>


void MixLdaU (int ns_size, double * new_ns_occ, double * ns_occ, std::unordered_map<std::string, InputKey *>& ControlMap, bool reset)
{
    RmgTimer RT0("Mix LdaU");

    if(Verify ("freeze_occupied", true, ControlMap)) return;

    /*Linear Mixing*/
    if (1 || Verify("charge_mixing_type","Linear", ControlMap) || ct.charge_pulay_order == 1 || ((ct.scf_steps < ct.davidson_premg) && (ct.md_steps == 0) && (ct.runflag != RESTART )) || (ct.xc_is_hybrid && Functional::is_exx_active()))
    {
        if(reset) return;
        RmgTimer RT1("Mix LdaU: Linear");
        for(int i = 0; i < ns_size; i++)
        {
            ns_occ[i] = ct.ldau_mix * new_ns_occ[i] + (1.0-ct.ldau_mix) * ns_occ[i]; 
        }
    }
    else 
    {
        RmgTimer RT1("Mix LdaU: Pulay");
        if(!Pulay_ldau)
        {
            Pulay_ldau = new PulayMixing(ns_size, ct.ldau_pulay_order, ct.ldau_pulay_refresh,
                    ct.ldau_mix, ct.ldau_pulay_scale, pct.grid_comm);
            Pulay_ldau->SetGspace(false, false, ct.drho_q0);

        }

        if(reset) {
            Pulay_ldau->Refresh();
            return;
        }

        double mone = -1.0;
        int ione = 1;
        daxpy(&ns_size, &mone, ns_occ, &ione, new_ns_occ, &ione);

        Pulay_ldau->Mixing(ns_occ, new_ns_occ);


    }
}

