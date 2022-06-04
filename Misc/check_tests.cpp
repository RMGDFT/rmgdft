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


#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h> 
#include <unistd.h>
#include <unordered_map>
#include <csignal>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "RmgThread.h"
#include "rmgthreads.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "Exxbase.h"
#include "Neb.h"
#include "Wannier.h"

void check_tests(void)
{
    FILE *elog = NULL, *blog = NULL;
    if(pct.imgpe == 0)
    {

        if(NULL == (elog = fopen ("test_energy.log", "w+")))
        {
            std::cout << "Unable to write test results to test_energy.log." << std::endl;
            return;
        }
        if(NULL == (blog = fopen ("test_bond_length.log", "w+")))
        {
            std::cout << "Unable to write test results to test_bond_length.log." << std::endl;
            return;
        }

        // Total energy test
        if(!std::isnan(ct.test_energy))
        {
            double eps = fabs(ct.TOTAL - ct.test_energy);
            fprintf(elog, "reference total energy: %15.8f Ha\n", ct.test_energy);
            fprintf(elog, "current   total energy: %15.8f Ha\n", ct.TOTAL);
            fprintf(elog, "deviation             : %15.8f Ha\n", eps);
            fprintf(elog, "tolerance             : %15.8f Ha\n", ct.test_energy_tolerance);
            if(eps < ct.test_energy_tolerance)
                fprintf(elog, "test status: pass\n");
            else
                fprintf(elog, "test status: fail\n");
        }

        // Bond length test. Always performed between the first and second atoms
        if(!std::isnan(ct.test_bond_length) && (Atoms.size() > 1))
        {
            double bond = 10000.0;
            for(int ix = -1; ix<= 1; ix++)
            {
                for(int iy = -1; iy<= 1; iy++)
                {
                    for(int iz = -1; iz<= 1; iz++)
                    {
                        double x = Atoms[0].crds[0] - Atoms[1].crds[0] + ix * Rmg_L.a0[0] + iy * Rmg_L.a1[0] + iz * Rmg_L.a2[0];
                        double y = Atoms[0].crds[1] - Atoms[1].crds[1] + ix * Rmg_L.a0[1] + iy * Rmg_L.a1[1] + iz * Rmg_L.a2[1];
                        double z = Atoms[0].crds[2] - Atoms[1].crds[2] + ix * Rmg_L.a0[2] + iy * Rmg_L.a1[2] + iz * Rmg_L.a2[2];
                        double r = sqrt(x*x + y*y + z*z);
                        if(r > 1.0e-5) bond = std::min(bond, r);
                    }
                }
            }

            double eps = fabs(ct.test_bond_length - bond);
            fprintf(blog, "reference bond length : %15.8f Bohr\n", ct.test_bond_length);
            fprintf(blog, "current   bond length : %15.8f Bohr\n", bond);
            fprintf(blog, "deviation             : %15.8f Bohr\n", eps);
            fprintf(blog, "tolerance             : %15.8f Bohr\n", ct.test_bond_length_tolerance);
            if(eps < ct.test_bond_length_tolerance)
                fprintf(blog, "test status: pass\n");
            else
                fprintf(blog, "test status: fail\n");
        }

        fclose(blog);
        fclose(elog);
    }
}

