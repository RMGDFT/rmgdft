/************************** SVN Revision Information **************************
 **    $Id: read_data.c 2426 2014-08-22 14:14:37Z ebriggs $    **
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
 *   void read_data(char *name, double *vh, double *rho, double *vxc, STATE *states)
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
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of doubles (doubles) but writes floats */
static void read_float (int fhand, double * rp, int count);
static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

template void ReadData(char *, double *, double *, double *, Kpoint<double> **);
template void ReadData(char *, double *, double *, double *, Kpoint<std::complex<double> > **);

/* Reads the hartree potential, the wavefunctions, the */
/* compensating charges and various other things from a file. */
template <typename KpointType>
void ReadData (char *name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr)
{
    char newname[MAX_PATH + 200];
    int grid[3];
    int fine[3];
    int pe[3];
    int npe;
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;

    /* wait until everybody gets here */
    /* my_barrier (); */
    MPI_Barrier(pct.img_comm);	

    /* Make the new output file name */
    rmg_printf("\nspin flag =%d\n", ct.spin_flag);
    if (ct.spin_flag)
    {
	if (pct.spinpe == 0)
	    sprintf(newname, "%s.up%d", name, pct.gridpe);
	else          /* if (pct.spinpe == 1)  */  
	    sprintf(newname, "%s.dw%d", name, pct.gridpe);
    }
    else
	sprintf (newname, "%s%d", name, pct.gridpe);


    int fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE);
    if (fhand < 0) {
        rmg_printf("Can't open data file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }



    /* read grid info */
    read_int (fhand, grid, 3);
    if (grid[0] != get_NX_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong NX_GRID");
    if (grid[1] != get_NY_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong NY_GRID");
    if (grid[2] != get_NZ_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong NZ_GRID");

    /* read grid processor topology */
    read_int (fhand, pe, 3);
    if (pe[0] != get_PE_X())
	rmg_error_handler (__FILE__, __LINE__,"Wrong PE_X");
    if (pe[1] != get_PE_Y())
	rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Y");
    if (pe[2] != get_PE_Z())
	rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Z");

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = (grid[0] * grid[1] * grid[2]) / npe;

    /* read fine grid info */
    read_int (fhand, fine, 3);
    if (fine[0] != get_FPX0_GRID() / get_PX0_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[1] != get_FPY0_GRID() / get_PY0_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[2] != get_FPZ0_GRID() / get_PZ0_GRID())
	rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];

    /* print out  */
    rmg_printf ("read_data: psi grid = %d %d %d\n", grid[0], grid[1], grid[2]);
    rmg_printf ("read_data: pe grid = %d %d %d\n", pe[0], pe[1], pe[2]);
    rmg_printf ("read_data: grid_size  = %d\n", grid_size);
    rmg_printf ("read_data: fine = %d %d %d\n", fine[0], fine[1], fine[2]);
    rmg_printf ("read_data: fgrid_size = %d\n", fgrid_size);


    /* read wavefunction info */
    read_int (fhand, &gamma, 1);
    if (gamma != ct.is_gamma)
	rmg_error_handler (__FILE__, __LINE__,"Wrong gamma data");


    read_int (fhand, &nk, 1);
    if (nk != ct.num_kpts && ct.forceflag != BAND_STRUCTURE)    /* bandstructure calculation */
	rmg_error_handler (__FILE__, __LINE__,"Wrong number of k points");

    rmg_printf ("read_data: gamma = %d\n", gamma);
    rmg_printf ("read_data: nk = %d\n", ct.num_kpts); 

    /* read number of states */  
    read_int (fhand, &ns, 1);
    if (ns != ct.num_states) {
	rmg_printf ("Wrong number of states: read %d from wave file, but ct.num_states is %d",ns, ct.num_states);
        rmg_error_handler (__FILE__, __LINE__,"Terminating.");
    }

    rmg_printf ("read_data: ns = %d\n", ns);


    /* read the hartree potential */
    read_double (fhand, vh, fgrid_size);
    rmg_printf ("read_data: read 'vh'\n");

    /* read density */
    read_double (fhand, rho, fgrid_size);
    rmg_printf ("read_data: read 'rho'\n");

    /* read Vxc */
    read_double (fhand, vxc, fgrid_size);
    rmg_printf ("read_data: read 'vxc'\n");


    /* read state occupations */
    {
	double *occ = new double[nk * ns];
        
	read_double (fhand, occ, (nk * ns));

	printf ("read_data: read 'occupations'\n"); 


	if (ct.forceflag != BAND_STRUCTURE)
	{
	    double occ_total = 0.0; 

	    for (ik = 0; ik < nk; ik++)
		for (is = 0; is < ns; is++)
		{
		    occ_total += ( Kptr[ik]->Kstates[is].occupation[0] = occ[ik * ns + is] );
		}



	    /* 
	       since we are using floats on the data file the precision
	       of the occupations is worse than 1e-10 required by the fill() routine
	       therefore we need to 'renormalize' the occupations so
	       that they add up to an integer
	       it's a hack I know, but whatever... not a biggie
	     */

	    {
		double iocc_total = (double) (int) (occ_total + 0.5);
		double fac = iocc_total / occ_total;

		    for (ik = 0; ik < nk; ik++)
			for (is = 0; is < ns; is++)
			{
			    Kptr[ik]->Kstates[is].occupation[0] *= fac;
			}
		/* end of normailization*/
	    }

	}             /* end if */

	delete [] occ;

    }           /* end of read occupations */




    /* read state eigenvalues, not needed really */
    {

	/* Read eigenvalue in pairwised case, while in polarized case, 
	 * it's the eigenvalue for proceesor's own spin  */ 
	for (ik = 0; ik < nk; ik++)
	    for (is = 0; is < ns; is++)
	    {
		read_double (fhand, &Kptr[ik]->Kstates[is].eig[0], 1);
	    }

	printf ("read_data: read 'eigenvalues'\n");

    }      /* end of read eigenvalues */



    /* read wavefunctions */
    {
	int wvfn_size = (gamma) ? grid_size : 2 * grid_size;

	for (ik = 0; ik < nk; ik++)
	{
	    for (is = 0; is < ns; is++)
	    {

		read_float (fhand, (double *)Kptr[ik]->Kstates[is].psi, wvfn_size);

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
/* The following routine accepts a buffer of doubles (doubles) but writes floats */

static void read_double (int fhand, double * rp, int count)
{

    ssize_t wanted = sizeof (double) * (ssize_t)count;
    ssize_t size = read (fhand, rp, wanted);
    if(size != wanted)
	rmg_error_handler (__FILE__, __LINE__,"error reading");

}
static void read_float (int fhand, double * rp, int count)
{

    int i;
    float *buf = new float[count];
    ssize_t wanted = sizeof (float) * (ssize_t)count;

    ssize_t size = read (fhand, buf, wanted);
    if(size != wanted)
	rmg_error_handler (__FILE__, __LINE__,"error reading");

    for (i = 0; i < count; i++)
	rp[i] = (double) buf[i];  /* floats take only 4 bytes instead of 8 bytes for double */

    delete [] buf;
}


static void read_int (int fhand, int *ip, int count)
{
    int size;

    size = count * sizeof (int);

    if (size != read (fhand, ip, size))
	rmg_error_handler (__FILE__, __LINE__,"error reading");

}


/******/
