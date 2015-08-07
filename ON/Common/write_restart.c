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
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "init_var.h"
#include "prototypes_on.h"
#include "transition.h"


void write_restart (char *name, double * vh, double *vxc, double *vh_old, double *vxc_old,  double * rho, STATE *states)
{
    char newname[MAX_PATH + 20];
    FILE *fhandle;
    int fhand;
    double time0, write_time;


    time0 = my_crtc ();
    
    /*If output file is specified as /dev/null or /dev/null/, skip writing */
    if ((!strcmp ("/dev/null", name)) || (!strcmp ("/dev/null/", name)) )
    {
	if (pct.imgpe == 0)
	    printf ("write_data: Output file given as /dev/null, no restart data written ...\n");
	return;
    }

    
    /*Only one processor will write restart file*/
    if (pct.imgpe == 0)
    {

	/*This opens restart file, creates a directory if needed */
	fhandle = open_restart_file (name);
	printf ("write_data: Restart file %s opened...\n", name);

	fclose (fhandle);

    } /*end if (pct.imgpe == 0) */

    /* All processors should wait until 0 is done to make sure that directories are created*/
    MPI_Barrier(pct.img_comm);

    if (ct.spin_flag)
    {   
	if (pct.spinpe==0)
	    sprintf (newname, "%s.up%d", name, pct.gridpe);
	else if(pct.spinpe==1) 
	    sprintf (newname, "%s.dw%d", name, pct.gridpe);

    }
    else
	sprintf (newname, "%s%d", name, pct.gridpe);


    write_data (name, vh, vxc, vh_old, vxc_old, rho, states);

    write_time = my_crtc () - time0;

    printf ("write_restart: writing took %.1f seconds \n", write_time);
    

    /* force change mode of output file */

}                               /* end write_data */

