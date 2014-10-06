/************************** SVN Revision Information **************************
 **    $Id: write_data.c 2426 2014-08-22 14:14:37Z ebriggs $    **
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
 *   void write_data(char *name, double *vh, double *rho, double *vxc, STATE *states)
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

static size_t totalsize;


/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of doubles (doubles) but writes floats */
static void write_float (int fh, double * rp, int count);
static void write_double (int fh, double * rp, int count);
static void write_int (int fh, int *ip, int count);


template void WriteData (int, double *, double *, double *, double *, Kpoint<double> **);
template void WriteData (int, double *, double *, double *, double *, Kpoint<std::complex<double> > **);


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
template <typename KpointType>
void WriteData (int fhand, double * vh, double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr)
{
    int fine[3];
    int grid[3];
    int pe[3];
    int npe;
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;
    double time0, write_time;

    time0 = my_crtc ();
    totalsize = 0;



    /* write grid info */
    grid[0] = get_NX_GRID();
    grid[1] = get_NY_GRID();
    grid[2] = get_NZ_GRID();
    write_int (fhand, grid, 3);

    /* write grid processor topology */
    pe[0] = get_PE_X();
    pe[1] = get_PE_Y();
    pe[2] = get_PE_Z();
    write_int (fhand, pe, 3);

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = (grid[0] * grid[1] * grid[2]) / npe;

    /* write fine grid info */
    fine[0] = get_FPX0_GRID() / get_PX0_GRID();
    fine[1] = get_FPY0_GRID() / get_PY0_GRID();
    fine[2] = get_FPZ0_GRID() / get_PZ0_GRID();
    write_int (fhand, fine, 3);
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];


    /* write wavefunction info */
    gamma = ct.is_gamma;
    nk = ct.num_kpts;
    write_int (fhand, &gamma, 1);
    write_int (fhand, &nk, 1); 

    ns = ct.num_states;
    write_int (fhand, &ns, 1); 


    /* write the hartree potential */
    write_double (fhand, vh, fgrid_size);

    /* write the total electronic density */
    write_double (fhand, rho, fgrid_size);

    /* write Vxc */
    write_double (fhand, vxc, fgrid_size);




    /* write the state occupations, in spin-polarized calculation, 
     * it's occupation for the processor's own spin */ 
    {
	for (ik = 0; ik < nk; ik++)
	    for (is = 0; is < ns; is++)
	    {
		write_double (fhand, &Kptr[ik]->Kstates[is].occupation[0], 1); 
	    }
    }
    

    /* write the state eigenvalues, while in spin-polarized case, 
     * it's eigenvalues of processor's own spin */
    {
	for (ik = 0; ik < nk; ik++)
	    for (is = 0; is < ns; is++)
	    {
		write_double (fhand, &Kptr[ik]->Kstates[is].eig[0], 1);
	    }

    }



    /* write wavefunctions */
    {
        int wvfn_size = (gamma) ? grid_size : 2 * grid_size;

        for (ik = 0; ik < nk; ik++)
        {
            for (is = 0; is < ns; is++)
            {
                write_float (fhand, (double *)Kptr[ik]->Kstates[is].psi, wvfn_size);
            }
        }
    }


    write_time = my_crtc () - time0;
    
    if (pct.imgpe == 0)
    {
        printf ("write_data: total size of each of the %d files = %.1f Mb\n", npe,
                ((double) totalsize) / (1024 * 1024));
        printf ("write_data: writing took %.1f seconds, writing speed %.3f Mbps \n", write_time,
                ((double) totalsize) / (1024 * 1024) / write_time);
    }




}                               /* end write_data */




/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of doubles (doubles) but writes floats */

static void write_double (int fh, double * rp, int count)
{
    int size;

    size = count * sizeof (double);
    if (size != write (fh, rp, size))
        rmg_error_handler (__FILE__,__LINE__,"error writing");

    totalsize += size;
}
static void write_float (int fh, double * rp, int count)
{
    float *buf = new float[count];
    int i, size;
    

    for (i = 0; i < count; i++)
        buf[i] = (float) rp[i]; /* 'float' uses only 4 bytes instead of 8 bytes for 'double' */

    size = count * sizeof (float);
    if (size != write (fh, buf, size))
        rmg_error_handler (__FILE__, __LINE__, "error writing");

    delete [] buf;
    totalsize += size;
}


static void write_int (int fh, int *ip, int count)
{
    int size;

    size = count * sizeof (int);
    if (size != write (fh, ip, size))
        rmg_error_handler (__FILE__, __LINE__, "error writing");

    totalsize += size;
}

/******/
