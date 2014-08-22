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
 *   void read_data(char *name, rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *vxc, STATE *states)
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
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"


/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of rmg_double_ts (doubles) but writes floats */
static void read_float (int fhand, rmg_double_t * rp, int count);
static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

/* Reads the hartree potential, the wavefunctions, the */
/* compensating charges and various other things from a file. */
void read_data (char *name, rmg_double_t * vh, rmg_double_t * rho, rmg_double_t * vxc, STATE * states)
{
    char newname[MAX_PATH + 200];
    int fhand;
    int grid[3];
    int fine[3];
    int pe[3];
    int npe;
    rmg_double_t a[9];
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;
    int na, ia;
    int i;
    rmg_double_t r[40];
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


    /* read grid info */
    read_int (fhand, grid, 3);
    if (grid[0] != get_NX_GRID())
	error_handler ("Wrong NX_GRID");
    if (grid[1] != get_NY_GRID())
	error_handler ("Wrong NY_GRID");
    if (grid[2] != get_NZ_GRID())
	error_handler ("Wrong NZ_GRID");

    /* read grid processor topology */
    read_int (fhand, pe, 3);
    if (pe[0] != get_PE_X())
	error_handler ("Wrong PE_X");
    if (pe[1] != get_PE_Y())
	error_handler ("Wrong PE_Y");
    if (pe[2] != get_PE_Z())
	error_handler ("Wrong PE_Z");

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = (grid[0] * grid[1] * grid[2]) / npe;

    /* read fine grid info */
    read_int (fhand, fine, 3);
    if (fine[0] != get_FPX0_GRID() / get_PX0_GRID())
	error_handler ("Wrong fine grid info");
    if (fine[1] != get_FPY0_GRID() / get_PY0_GRID())
	error_handler ("Wrong fine grid info");
    if (fine[2] != get_FPZ0_GRID() / get_PZ0_GRID())
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
    if (gamma != ct.is_gamma)
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
	rmg_double_t *occ;
	my_malloc (occ, nk * ns, rmg_double_t); 
	read_double (fhand, occ, (nk * ns));

	printf ("read_data: read 'occupations'\n"); 


	if (ct.forceflag != BAND_STRUCTURE)
	{
	    rmg_double_t occ_total = 0.0; 

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
		rmg_double_t iocc_total = (rmg_double_t) (int) (occ_total + 0.5);
		rmg_double_t fac = iocc_total / occ_total;

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


    close (fhand);


}                               /* end read_data */



/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of rmg_double_ts (doubles) but writes floats */

static void read_double (int fhand, double * rp, int count)
{

    int i, size;

    size = count * sizeof (double);

    if (size != read (fhand, rp, size))
	error_handler ("error reading");

}
static void read_float (int fhand, rmg_double_t * rp, int count)
{

    float *buf;
    int i, size;
    my_malloc (buf, count, float);

    size = count * sizeof (float);

    if (size != read (fhand, buf, size))
	error_handler ("error reading");

    for (i = 0; i < count; i++)
	rp[i] = (rmg_double_t) buf[i];  /* floats take only 4 bytes instead of 8 bytes for double */

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
