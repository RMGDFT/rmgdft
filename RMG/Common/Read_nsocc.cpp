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
#include <unistd.h>
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"
#include "ZfpCompress.h"
#include "RmgGemm.h"
#include "RmgException.h"

template void Read_nsocc(char *, Kpoint<double> *);
template void Read_nsocc(char *, Kpoint<std::complex<double> > *);

/* Reads the hartree potential, the wavefunctions, the */
/* compensating charges and various other things from a file. */
template <typename KpointType>
void Read_nsocc(char *name, Kpoint<KpointType> * kptr)
{
    char newname[MAX_PATH + 200];
    int pstride = kptr->ldaU->ldaU_m;
    size_t occ_size_bytes = ct.nspin * Atoms.size() * pstride * pstride * sizeof(std::complex<double>);

    if(pct.imgpe == 0)
    {
        sprintf (newname, "%s_spin%d_nsocc", name, pct.spinpe);
        int fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE);
        if (fhand < 0) {
            rmg_printf("Can't open data file %s", newname);
            rmg_error_handler(__FILE__, __LINE__, "Terminating.");
        }


        read(fhand, kptr->ldaU->ns_occ.data(), occ_size_bytes);
        close(fhand);
    }

    int occ_size = ct.nspin * Atoms.size() * pstride * pstride;
    MPI_Bcast(kptr->ldaU->ns_occ.data(), occ_size*2, MPI_DOUBLE, 0, pct.img_comm);

}



