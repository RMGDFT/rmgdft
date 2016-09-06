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

template void WriteRestart (char *, double *, double *, double *, double *, Kpoint<double> **);
template void WriteRestart (char *, double *, double *, double *, double *, Kpoint<std::complex<double> > **);

/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
template <typename KpointType>
void WriteRestart (char *name, double * vh, double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr)
{
    char newname[MAX_PATH + 20];
    int amode;
    FILE *fhandle;
    int fhand;
    int ia, ion;
    double time0, write_time;
    SPECIES *sp;
    ION *iptr;


    time0 = my_crtc ();
    
    /*If output file is specified as /dev/null or /dev/null/, skip writing */
    if ((!strcmp ("/dev/null", name)) || (!strcmp ("/dev/null/", name)) )
    {
	if (pct.imgpe == 0)
	    rmg_printf ("WriteRestart: Output file given as /dev/null, no restart data written ...\n");
	return;
    }

    
    /*Only one processor will write restart file*/
    if (pct.imgpe == 0)
    {

	/*This opens restart file, creates a directory if needed */
	fhandle = open_restart_file (name);
	rmg_printf ("WriteRestart: Restart file %s opened...\n", name);


        // Absolute coordinates in bohr
        fprintf(fhandle, "atomic_coordinate_type = \"Absolute\"\n");
        fprintf(fhandle, "crds_units = \"Bohr\"\n");

	/* write current ionic cartesian positions */
	fprintf(fhandle,"atoms = \"");

	for (ion = 0; ion < ct.num_ions; ion++)
	{

	    iptr = &ct.ions[ion];
	    sp = &ct.sp[iptr->species];

	    fprintf(fhandle,"\n %s %#15.12g %#15.12g %#15.12g %d", sp->atomic_symbol, iptr->crds[0], iptr->crds[1], iptr->crds[2], iptr->movable);
	}

	fprintf(fhandle,"\n\"\n");


	/* write current ionic velocities */
	fprintf(fhandle,"\nionic_velocities = \"");
	for (ion = 0; ion < ct.num_ions; ion++)
	{
	    iptr = &ct.ions[ion];
	    fprintf(fhandle, "\n %#15.12g %#15.12g %#15.12g ", iptr->velocity[0], iptr->velocity[1], iptr->velocity[2]);
	
	}
	fprintf(fhandle,"\n\"\n");

	/* write current ionic forces pointer array */
	fprintf(fhandle,"\nionic_force_pointer = \"%d %d %d %d\"\n", ct.fpt[0], ct.fpt[1], ct.fpt[2], ct.fpt[3]);

	/* write current ionic forces */
	fprintf(fhandle,"\nionic_forces = \"");
        for(int ic = 0;ic < 4;ic++) {
            for (ion = 0; ion < ct.num_ions; ion++)
            {
                iptr = &ct.ions[ion];
                fprintf(fhandle, "\n %#15.12g %#15.12g %#15.12g ", iptr->force[ic][0], iptr->force[ic][1], iptr->force[ic][2]);
            }
        }
	fprintf(fhandle,"\n\"\n");


	fprintf(fhandle,"\nnose_positions = \"%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\"\n", 
		ct.nose.xx[0], ct.nose.xx[1], ct.nose.xx[2], ct.nose.xx[3], ct.nose.xx[4], 
		ct.nose.xx[5], ct.nose.xx[6], ct.nose.xx[7], ct.nose.xx[8], ct.nose.xx[9]);

	fprintf(fhandle,"\nnose_velocities = \"%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\"\n", 
		ct.nose.xv[0], ct.nose.xv[1], ct.nose.xv[2], ct.nose.xv[3], ct.nose.xv[4], 
		ct.nose.xv[5], ct.nose.xv[6], ct.nose.xv[7], ct.nose.xv[8], ct.nose.xx[9]);

	fprintf(fhandle,"\nnose_masses = \"%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\"\n", 
		ct.nose.xq[0], ct.nose.xq[1], ct.nose.xq[2], ct.nose.xq[3], ct.nose.xq[4], 
		ct.nose.xq[5], ct.nose.xq[6], ct.nose.xq[7], ct.nose.xq[8], ct.nose.xq[9]);

	fprintf(fhandle,"\nnose_forces = \"");
	for (ia=0; ia<4; ia++)
	    fprintf(fhandle,"\n%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g", 
		    ct.nose.xf[ia][0], ct.nose.xf[ia][1], ct.nose.xf[ia][2], ct.nose.xf[ia][3], ct.nose.xf[ia][4], 
		    ct.nose.xf[ia][5], ct.nose.xf[ia][6], ct.nose.xf[ia][7], ct.nose.xf[ia][8], ct.nose.xf[ia][9]);
	fprintf(fhandle,"\n\"\n");

	fprintf(fhandle,"\nionic_time_step = \"%.12g\"", ct.iondt);
	fprintf(fhandle,"\ndynamic_time_counter = \"%d\"", ct.relax_steps_counter);
	fprintf(fhandle,"\nkpoint_distribution = \"%d\"", pct.pe_kpoint);
//	fprintf(fhandle,"\n");


	/* done with writing */
	fclose (fhandle);

    } /*end if (pct.imgpe == 0) */

    /* All processors should wait until 0 is done to make sure that directories are created*/
    MPI_Barrier(pct.img_comm);

    sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", name, pct.spinpe, pct.kstart, pct.gridpe);

    amode = S_IREAD | S_IWRITE;
    fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);
    if (fhand < 0) {
        rmg_printf("Can't open restart file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }

    WriteData (fhand, vh, rho, rho_oppo, vxc, Kptr);
    close (fhand);

    write_time = my_crtc () - time0;

    rmg_printf ("WriteRestart: writing took %.1f seconds \n", write_time);
    

    /* force change mode of output file */
    chmod (newname, amode);

}                               /* end write_data */

