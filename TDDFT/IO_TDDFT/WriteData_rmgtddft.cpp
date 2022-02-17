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
    #include <libgen.h>
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


void WriteData_rmgtddft (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *Smatrix, 
        double *Cmatrix, double *Hmatrix_m1, double *Hmatrix_0, int tot_steps, int n2)
{
    int fhand;
    int fgrid_size;
    char newname[MAX_PATH], tmpname[MAX_PATH];


    int amode;
    sprintf (newname, "%s%d", filename, pct.gridpe);

    amode = S_IREAD | S_IWRITE;
    fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);

   if (fhand < 0) {

       strcpy (tmpname, newname);


#if !(defined(_WIN32) || defined(_WIN64))
       if (!(mkdir (dirname (tmpname), S_IRWXU)))
#else
           char dirname[_MAX_DIR];
       _splitpath(tmpname, NULL, dirname, NULL, NULL);
       if (!_mkdir(dirname));
#endif
        rmg_printf ("\n Creating directory %s succesfull\n\n", tmpname);
    else
        rmg_printf ("\n Creating directory %s FAILED\n\n", tmpname);


    fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);

    if (fhand < 0) {
        rmg_printf("Can't open restart file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }


   }


   fgrid_size = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
   write (fhand, vh, fgrid_size * sizeof(double));
   write (fhand, vxc, fgrid_size * sizeof(double));
   write (fhand, vh_corr, fgrid_size * sizeof(double));

   write (fhand, Pn0, 2* n2 * sizeof(double));
   write (fhand, Hmatrix, n2 * sizeof(double));
   write (fhand, Smatrix, n2 * sizeof(double));
   write (fhand, Cmatrix, n2 * sizeof(double));
   write (fhand, Hmatrix_m1, n2 * sizeof(double));
   write (fhand, Hmatrix_0, n2 * sizeof(double));
   write (fhand, &tot_steps, sizeof(int));
   close(fhand);

}                               /* end write_data */
