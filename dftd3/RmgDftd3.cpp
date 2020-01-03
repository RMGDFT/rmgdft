
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

#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>
#include <algorithm>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgTimer.h"
#include "RmgSumAll.h"
#include "transition.h"
#include "MapElements.h"
#include "Functional.h"

#include "rmg_mangling.h"

#define dftd3_driver       RMG_FC_GLOBAL(dftd3_driver, DFTD3_DRIVER)

extern "C" void dftd3_driver_(int *nAtoms, int *nSpecies, double *coords, int *atnums, double  *latt, const char *functype, 
    int *d3_version, double *edisp, double *grads, double *stress, int *func_len);
void RmgDftd3(double *edisp, double *grads, double *stress, int d3_version) 
{
    
    int nAtoms = (int)Atoms.size();
    int nSpecies = ct.num_species;
    double *coords = new double[3 * Atoms.size()];
    int *atnums = new int[Atoms.size()];
   
    const std::string  functype = Functional::get_dft_name_rmg();
    int func_len  = functype.size();
    double latt[9];

    for(int ion = 0; ion < nAtoms; ion++)
    {
        coords[ion * 3 + 0] = Atoms[ion].crds[0];
        coords[ion * 3 + 1] = Atoms[ion].crds[1];
        coords[ion * 3 + 2] = Atoms[ion].crds[2];
        atnums[ion] = GetAtomicNumber(Atoms[ion].symbol);
    }

    for(int j = 0; j < 3; j++)
    {
        latt[0*3 + j] = Rmg_L.a0[j];
        latt[1*3 + j] = Rmg_L.a1[j];
        latt[2*3 + j] = Rmg_L.a2[j];
    }

    dftd3_driver(&nAtoms, &nSpecies, coords, atnums, latt, functype.c_str(), &d3_version, edisp, grads, stress, &func_len);
    delete [] coords;
    delete [] atnums;
}


