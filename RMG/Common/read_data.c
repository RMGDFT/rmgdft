/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/read_data.c *****
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
 *   void read_data(char *name, REAL *vh, REAL *rho, REAL *vxc, STATE *states)
 *   when ct.runflag == 1,
 *   Reads the hartree potential, the wavefunctions
 *   and various other things from a file which is created by the 
 *   previous run.          
 * INPUTS
 *   name:  file name
 * OUTPUT
 *   vh: Hartree potential
 *   rho:  total valence charge density
 *   vxc:  exchange correlation potential
 *   states: point to orbital structure
 * PARENTS
 *   init.c
 * CHILDREN
 *   scatter_psi.c
 * SEE ALSO
 *   write_data.c
 * SOURCE
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "main.h"


/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of REALs (doubles) but writes floats */
static void read_float (int fhand, REAL * rp, int count);
static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

/* Reads the hartree potential, the wavefunctions, the */
/* compensating charges and various other things from a file. */
void read_data (char *name, REAL * vh, REAL * rho, REAL * vxc, STATE * states)
{
    char newname[MAX_PATH + 200];
    int fhand;
    int grid[3];
    int fine[3];
    int pe[3];
    int npe;
    REAL a[9];
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;
    int na, ia;
    int i;
    REAL r[40];
    int tmp_int[4];

    /* wait until everybody gets here */
    /* my_barrier (); */
    MPI_Barrier(pct.img_comm);	

    sprintf(newname, "%s.restart", name);
    read_control(newname);

    /* Make the new output file name */
    printf("\nspin flag =%d\n", ct.spin_flag);
    if (ct.spin_flag)
    {
	if (pct.spinpe == 0)
	    sprintf(newname, "%s.up%d", name, pct.gridpe);
	else          /* if (pct.spinpe == 1)  */  
	    sprintf(newname, "%s.dw%d", name, pct.gridpe);
    }
    else
	sprintf (newname, "%s%d", name, pct.gridpe);


    my_open (fhand, newname, O_RDWR, S_IREAD | S_IWRITE);

    /* read lattice info */
    read_double (fhand, a, 9);

    /* read grid info */
    read_int (fhand, grid, 3);
    if (grid[0] != NX_GRID)
	error_handler ("Wrong NX_GRID");
    if (grid[1] != NY_GRID)
	error_handler ("Wrong NY_GRID");
    if (grid[2] != NZ_GRID)
	error_handler ("Wrong NZ_GRID");

    /* read grid processor topology */
    read_int (fhand, pe, 3);
    if (pe[0] != PE_X)
	error_handler ("Wrong PE_X");
    if (pe[1] != PE_Y)
	error_handler ("Wrong PE_Y");
    if (pe[2] != PE_Z)
	error_handler ("Wrong PE_Z");

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = (grid[0] * grid[1] * grid[2]) / npe;

    /* read fine grid info */
    read_int (fhand, fine, 3);
    if (fine[0] != FPX0_GRID / PX0_GRID)
	error_handler ("Wrong fine grid info");
    if (fine[1] != FPY0_GRID / PY0_GRID)
	error_handler ("Wrong fine grid info");
    if (fine[2] != FPZ0_GRID / PZ0_GRID)
	error_handler ("Wrong fine grid info");
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];

    /* print out  */
    printf ("read_data: psi grid = %d %d %d\n", grid[0], grid[1], grid[2]);
    printf ("read_data: pe grid = %d %d %d\n", pe[0], pe[1], pe[2]);
    printf ("read_data: grid_size  = %d\n", grid_size);
    printf ("read_data: fine = %d %d %d\n", fine[0], fine[1], fine[2]);
    printf ("read_data: fgrid_size = %d\n", fgrid_size);




    /* read wavefunction info */
    read_int (fhand, &gamma, 1);
    if (gamma != GAMMA_PT)
	error_handler ("Wrong gamma data");


    read_int (fhand, &nk, 1);
    if (nk != ct.num_kpts && ct.forceflag != BAND_STRUCTURE)    /* bandstructure calculation */
	error_handler ("Wrong number of k points");

    printf ("read_data: gamma = %d\n", gamma);
    printf ("read_data: nk = %d\n", ct.num_kpts); 

    /* read number of states */  
    read_int (fhand, &ns, 1);
    if (ns != ct.num_states)
	error_handler ("Wrong number of states: read %d from wave file, but ct.num_states is %d",ns, ct.num_states);

    printf ("read_data: ns = %d\n", ns);


    /* read the hartree potential */
    read_double (fhand, vh, fgrid_size);
    printf ("read_data: read 'vh'\n");

    /* read density */
    read_double (fhand, rho, fgrid_size);
    printf ("read_data: read 'rho'\n");

    /* read Vxc */
    read_double (fhand, vxc, fgrid_size);
    printf ("read_data: read 'vxc'\n");


    /* read state occupations */
    {
	STATE *sp;
	REAL *occ;
	my_malloc (occ, nk * ns, REAL); 
	read_double (fhand, occ, (nk * ns));

	printf ("read_data: read 'occupations'\n"); 


	if (ct.forceflag != BAND_STRUCTURE)
	{
	    REAL occ_total = 0.0; 

	    sp = states;
	    for (ik = 0; ik < nk; ik++)
		for (is = 0; is < ns; is++)
		{
		    occ_total += ( sp->occupation[0] = occ[ik * ns + is] );
		    sp++;
		}



	    /* 
	       since we are using floats on the data file the precision
	       of the occupations is worse than 1e-10 required by the fill() routine
	       therefore we need to 'renormalize' the occupations so
	       that they add up to an integer
	       it's a hack I know, but whatever... not a biggie
	     */

	    {
		REAL iocc_total = (REAL) (int) (occ_total + 0.5);
		REAL fac = iocc_total / occ_total;

		    sp = states;
		    for (ik = 0; ik < nk; ik++)
			for (is = 0; is < ns; is++)
			{
			    sp->occupation[0] *= fac;
			    sp++;
			}
		/* end of normailization*/
	    }

	}             /* end if */

	my_free (occ);

    }           /* end of read occupations */




    /* read state eigenvalues, not needed really */
    {

	STATE *sp;

	/* Read eigenvalue in pairwised case, while in polarized case, 
	 * it's the eigenvalue for proceesor's own spin  */ 
	sp = states;
	for (ik = 0; ik < nk; ik++)
	    for (is = 0; is < ns; is++)
	    {
		read_double (fhand, &sp->eig[0], 1);
		sp++;
	    }

	printf ("read_data: read 'eigenvalues'\n");

    }      /* end of read eigenvalues */



    /* read wavefunctions */
    {
	int wvfn_size = (gamma) ? grid_size : 2 * grid_size;
	STATE *sp;

	sp = states;
	for (ik = 0; ik < nk; ik++)
	{
	    for (is = 0; is < ns; is++)
	    {

		read_double (fhand, sp->psiR, wvfn_size);
		sp++;

	    }
	    /*  for calculating band structures, 
		only read-in one wave function (?) */
	    if (ct.forceflag == BAND_STRUCTURE)
		return;
	}

	printf ("read_data: read 'wfns'\n");

    }


    /* read number of ions */
    read_int (fhand, &na, 1);
    /*if (na != ct.num_ions)
	error_handler ("Wrong number of ions");*/


    /* read current ionic cartesian positions */
    {
	for (ia = 0; ia < na; ia++)
	{
	    read_double (fhand, r, 3);

	    /*for (i = 0; i < 3; i++)
		ct.ions[ia].crds[i] = r[i];*/

	}


	/* read current ionic crystal positions */
	for (ia = 0; ia < na; ia++)
	{
	    read_double (fhand, r, 3);

	    /*for (i = 0; i < 3; i++)
		ct.ions[ia].xtal[i] =
		    (r[i] < 0) ? (r[i] + 1.0) : ((r[i] > 1.0) ? (r[i] - 1.0) : r[i]);*/

	    /*to_cartesian (ct.ions[ia].xtal, ct.ions[ia].crds);*/

	}

	/* Overwrite the initial positions with current positions, 
	 * may be useful for constraint dynamics  */

	/* read original ionic cartesian positions */
	for (ia = 0; ia < na; ia++)
	{
	    //read_double (fhand, &ct.ions[ia].icrds[0], 3);
	    read_double (fhand, r, 3);

	}


	/* read original ionic crystal positions */
	for (ia = 0; ia < na; ia++)
	{
	    //read_double (fhand, &ct.ions[ia].ixtal[0], 3);
	    read_double (fhand, r, 3);

	}


	/* read ionic velocities */
	for (ia = 0; ia < na; ia++)
	{
	    //read_double (fhand, &ct.ions[ia].velocity[0], 3);
	    read_double (fhand, r, 3);
	}

	/* read forces pointer */
	//read_int (fhand, &ct.fpt[0], 4);
	read_int (fhand, tmp_int, 4);

	/* read ionic forces */
	/*for (ia = 0; ia < na; ia++)
	    read_double (fhand, &ct.ions[ia].force[0][0], 3 * 4);*/
	for (ia = 0; ia < 4*na; ia++)
	    read_double (fhand, r, 3);


	/* read Nose positions,velocities, masses and forces from the file */
	read_double (fhand, r, 10);
	read_double (fhand, r, 10);
	read_double (fhand, r, 10);
	read_double (fhand, r, 4 * 10);

	/* read ionic timestep*/
	read_double (fhand, r, 1);

    }

    close (fhand);


}                               /* end read_data */



/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of REALs (doubles) but writes floats */

static void read_double (int fhand, double * rp, int count)
{

    int i, size;

    size = count * sizeof (double);

    if (size != read (fhand, rp, size))
	error_handler ("error reading");

}
static void read_float (int fhand, REAL * rp, int count)
{

    float *buf;
    int i, size;
    my_malloc (buf, count, float);

    size = count * sizeof (float);

    if (size != read (fhand, buf, size))
	error_handler ("error reading");

    for (i = 0; i < count; i++)
	rp[i] = (REAL) buf[i];  /* floats take only 4 bytes instead of 8 bytes for double */

    my_free (buf);
}


static void read_int (int fhand, int *ip, int count)
{
    int size;

    size = count * sizeof (int);

    if (size != read (fhand, ip, size))
	error_handler ("error reading");

}


/******/
