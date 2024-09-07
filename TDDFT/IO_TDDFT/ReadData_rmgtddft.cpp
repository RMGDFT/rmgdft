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

void ReadData_rmgtddft (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix,  
        double *Hmatrix_m1, double *Hmatrix_0, int *tot_steps, int n2,
        std::vector<double> &Eterms, double *Hcore_tddft, int numst)
{
    int fhand, fgrid_size;
    char newname[MAX_PATH];


    int amode;
    sprintf (newname, "%s%d", filename, pct.gridpe);

    amode = S_IREAD | S_IWRITE;
    fhand = open(newname, O_RDWR, amode);

   if (fhand < 0) {
        rmg_printf("Can't open restart file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }


    fgrid_size = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
    read (fhand, vh, fgrid_size * sizeof(double));
    read (fhand, vxc, fgrid_size * sizeof(double));
    read (fhand, vh_corr, fgrid_size * sizeof(double));
    double *PnI = Pn0 + n2;
    read (fhand, Pn0, n2 * sizeof(double));
    read (fhand, PnI, n2 * sizeof(double));
    read (fhand, Hmatrix, n2 * sizeof(double));
    read (fhand, Hmatrix_m1, n2 * sizeof(double));
    read (fhand, Hmatrix_0, n2 * sizeof(double));
    read (fhand, tot_steps, sizeof(int));
    read (fhand, Eterms.data(), Eterms.size() * sizeof(double) );
    size_t size = read (fhand, Hcore_tddft, numst * numst * sizeof(double));


    if(size != numst * numst * sizeof(double)) 
        rmg_error_handler(__FILE__, __LINE__, "endof file in ReadData_rmgtddft ");

    close(fhand);

}                               /* end write_data */
