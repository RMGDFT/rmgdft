/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/write_data.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void write_data(char *name, rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *vxc, STATE *states)
 *   Writes the hartree potential, the wavefunctions, the 
 *   charge density and various other things to a file.
 *   This file is useful for re-run (ct.runflag =1)
 * INPUTS
 *   name: file name
 *   vh:  Hartree potential
 *   rho: total valence charge density
 *   vxc: exchange correlation potential
 *   states: points to orbital structure
 * OUTPUT
 *   write to a file 
 * PARENTS
 *   main.c
 * CHILDREN
 *   gather_psi.c 
 * SEE ALSO
 *   read_data.c
 * SOURCE
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "main.h"


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void write_restart (char *name, rmg_double_t * vh, rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * vxc, STATE * states)
{
    char newname[MAX_PATH + 20];
    int amode;
    FILE *fhandle;
    int fhand;
    int ia, ion;
    rmg_double_t time0, write_time;
    SPECIES *sp;
    ION *iptr;


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
	for (ion = 0; ion < ct.num_ions; ion++)
	{
	    iptr = &ct.ions[ion];
	    fprintf(fhandle, "\n %#15.12g %#15.12g %#15.12g ", iptr->force[0][0], iptr->force[0][1], iptr->force[0][2]);
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
//	fprintf(fhandle,"\n");


	/* done with writing */
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


    my_open (fhand, newname, O_CREAT | O_TRUNC | O_RDWR, amode);
    write_data (fhand, vh, rho, rho_oppo, vxc, states);
    close (fhand);

    write_time = my_crtc () - time0;

     printf ("write_restart: writing took %.1f seconds \n", write_time);
    

    /* force change mode of output file */
    amode = S_IREAD | S_IWRITE;
    chmod (newname, amode);

}                               /* end write_data */

/******/
