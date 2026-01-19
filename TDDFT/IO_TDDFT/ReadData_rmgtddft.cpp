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
#include "State.h"
#include "Kpoint.h"
#include "transition.h"

void ReadData_rmgtddft (const char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix,  
        double *Hmatrix_m1, double *Hmatrix_0, int *tot_steps, int n2, int n2_C,
        std::vector<double> &Eterms, double *Hcore_tddft, int numst)
{
    int fhand, fgrid_size;

    int amode;
    amode = S_IREAD | S_IWRITE;
    fhand = open(filename, O_RDWR, amode);

   if (fhand < 0) {
        rmg::printlog("Can't open restart file %s", filename);
        rmg::error("Terminating.");
    }


    int n_rho = ct.noncoll_factor * ct.noncoll_factor;

    fgrid_size = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
    rmg::readfile (fhand, vh, fgrid_size * sizeof(double));
    rmg::readfile (fhand, vxc, n_rho*fgrid_size * sizeof(double));
    rmg::readfile (fhand, vh_corr, fgrid_size * sizeof(double));
    rmg::readfile (fhand, Pn0, 2*n2 * sizeof(double));
    rmg::readfile (fhand, Hmatrix, n2_C * sizeof(double));
    rmg::readfile (fhand, Hmatrix_m1, n2_C * sizeof(double));
    rmg::readfile (fhand, Hmatrix_0, n2_C * sizeof(double));
    rmg::readfile (fhand, tot_steps, sizeof(int));
    rmg::readfile (fhand, Eterms.data(), Eterms.size() * sizeof(double) );
    rmg::readfile (fhand, Hcore_tddft, numst * numst * sizeof(double));
    close(fhand);

}                               /* end read_data */
