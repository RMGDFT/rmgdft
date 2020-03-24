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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <unistd.h>
#else
    #include <io.h>
#endif
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "LocalObject.h"
#include "prototypes_on.h"
#include "blas.h"



/* 

  This routine writes the charge density, hartree and exchange correlation potentials and
  the electronic orbitals to serial files where each global grid object represents one file.

*/

template void WriteCij (std::string&, double *Cij);
//template void WriteWavefunctions (std::string&, LocalObjectt<std::complex<double> > &LO, std::complex<double> &Cij);

template <typename KpointType>
void WriteCij (std::string& name, KpointType *Cij_dist)
{
    KpointType *Cij_global = new KpointType[ct.num_states * ct.num_states]();

    mat_dist_to_global(Cij_dist, pct.desca, Cij_global);

    std::string wfname = name + "_spin" + std::to_string(pct.spinpe) + "_Cij" ;

    if(pct.gridpe == 0)
    {
        int fhand = open(wfname.c_str(), O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE);
        if (fhand < 0) {
            rmg_printf("Can't open restart file %s", wfname.c_str());
            rmg_error_handler(__FILE__, __LINE__, "Terminating.");
        }
        size_t size = ct.num_states * ct.num_states * sizeof(double);
        write (fhand, Cij_global, size);
        close(fhand);
        fflush(NULL);
    }


    delete [] Cij_global;

} // WriteSerialData

