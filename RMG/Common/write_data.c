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
 *   void write_data(char *name, REAL *vh, REAL *rho, REAL *vxc, STATE *states)
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

static size_t totalsize;


/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of REALs (doubles) but writes floats */
static void write_float (int fh, REAL * rp, int count);
static void write_double (int fh, double * rp, int count);
static void write_int (int fh, int *ip, int count);




/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void write_data (int fhand, REAL * vh, REAL * rho, REAL * rho_oppo, REAL * vxc, STATE * states)
{
    char newname[MAX_PATH + 20];
    int amode;
    int fine[3];
    int grid[3];
    int pe[3];
    int npe;
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;
    int ia, idx, nspin = (ct.spin_flag + 1);
    REAL time0, write_time;

    time0 = my_crtc ();
    totalsize = 0;


    /* write lattice information */
    write_double (fhand, ct.a0, 3);
    write_double (fhand, ct.a1, 3);
    write_double (fhand, ct.a2, 3);

    /* write grid info */
    grid[0] = NX_GRID;
    grid[1] = NY_GRID;
    grid[2] = NZ_GRID;
    write_int (fhand, grid, 3);

    /* write grid processor topology */
    pe[0] = PE_X;
    pe[1] = PE_Y;
    pe[2] = PE_Z;
    write_int (fhand, pe, 3);

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = (grid[0] * grid[1] * grid[2]) / npe;

    /* write fine grid info */
    fine[0] = FPX0_GRID / PX0_GRID;
    fine[1] = FPY0_GRID / PY0_GRID;
    fine[2] = FPZ0_GRID / PZ0_GRID;
    write_int (fhand, fine, 3);
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];


    /* write wavefunction info */
    gamma = GAMMA_PT;
    nk = ct.num_kpts;
    write_int (fhand, &gamma, 1);
    write_int (fhand, &nk, 1); 

    ns = ct.num_states;
    write_int (fhand, &ns, 1); 


    /* write the hartree potential */
    write_double (fhand, vh, fgrid_size);

    /* write the total electronic density */
    write_double (fhand, rho, fgrid_size);

    if (ct.spin_flag)
    {
    	/* write the electronic density for the opposite spin */
   	 write_double (fhand, rho_oppo, fgrid_size); 
    }

    /* write Vxc */
    write_double (fhand, vxc, fgrid_size);




    /* write the state occupations, in spin-polarized calculation, 
     * it's occupation for the processor's own spin */ 
    {
        STATE *sp;
	for (idx = 0; idx < nspin; idx++)
	{
            sp = states;
            for (ik = 0; ik < nk; ik++)
            	for (is = 0; is < ns; is++)
            	{
                	write_double (fhand, &sp->occupation[idx], 1); 
                	sp++;
            	}
	}
    }
    

    /* write the state eigenvalues, while in spin-polarized case, 
     * it's eigenvalues of processor's own spin */
    {
        STATE *sp;
	for (idx = 0; idx < nspin; idx++)
	{
        	sp = states;
                for (ik = 0; ik < nk; ik++)
            		for (is = 0; is < ns; is++)
            		{
                		write_double (fhand, &sp->eig[idx], 1);
                		sp++;
            		}
	}

    }



    /* write wavefunctions */
    {
        int wvfn_size = (gamma) ? grid_size : 2 * grid_size;
        STATE *sp;


        sp = states;
        for (ik = 0; ik < nk; ik++)
        {
            for (is = 0; is < ns; is++)
            {
                write_double (fhand, sp->psiR, wvfn_size);
                sp++;
            }
        }
    }





    /* write number of ions */
    write_int (fhand, &ct.num_ions, 1);

    /* write current ionic cartesian positions */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, ct.ions[ia].crds, 3);


    /* write current ionic crystal positions */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, ct.ions[ia].xtal, 3);

    /* write original ionic cartesian positions */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, ct.ions[ia].icrds, 3);

    /* write original ionic crystal positions */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, ct.ions[ia].ixtal, 3);


    /* write current ionic velocities */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, ct.ions[ia].velocity, 3);

    /* write current ionic forces pointer array */
    write_int (fhand, ct.fpt, 4);

    /* write current ionic forces */
    for (ia = 0; ia < ct.num_ions; ia++)
        write_double (fhand, &ct.ions[ia].force[0][0], 3 * 4);

    /* write Nose positions,velocities, masses and forces */
    write_double (fhand, &ct.nose.xx[0], 10);
    write_double (fhand, &ct.nose.xv[0], 10);
    write_double (fhand, &ct.nose.xq[0], 10);
    write_double (fhand, &ct.nose.xf[0][0], 4 * 10);

    /*Write ionic timestep*/
    write_double (fhand, &ct.iondt, 1);
    
    write_time = my_crtc () - time0;
    
    if (pct.imgpe == 0)
    {
        /*Things like grids and number of precessors should be written in the beginning of the file */
        /*printf( "write_data: psi grid = %d %d %d\n", ct.psi_nxgrid, ct.psi_nygrid, ct.psi_nzgrid);
           printf( "write_data: pe grid = %d %d %d\n", PE_X, PE_Y, PE_Z);
           printf( "write_data: grid_size = %d\n", grid_size);
           printf( "write_data: gamma = %d\n", gamma);
           printf( "write_data: nk = %d\n", nk);
           printf( "write_data: ns = %d\n", ns); */
        printf ("write_data: total size of each of the %d files = %.1f Mb\n", npe,
                ((REAL) totalsize) / (1024 * 1024));
        printf ("write_data: writing took %.1f seconds, writing speed %.3f Mbps \n", write_time,
                ((REAL) totalsize) / (1024 * 1024) / write_time);
    }




}                               /* end write_data */




/* To save disk space 'floats' are written instead of 'doubles'. */
/* The following routine accepts a buffer of REALs (doubles) but writes floats */

static void write_double (int fh, double * rp, int count)
{
    int size;

    size = count * sizeof (double);
    if (size != write (fh, rp, size))
        error_handler ("error writing");

    totalsize += size;
}
static void write_float (int fh, REAL * rp, int count)
{
    float *buf;
    int i, size;
    my_malloc (buf, count, float);

    for (i = 0; i < count; i++)
        buf[i] = (float) rp[i]; /* 'float' uses only 4 bytes instead of 8 bytes for 'double' */

    size = count * sizeof (float);
    if (size != write (fh, buf, size))
        error_handler ("error writing");

    my_free (buf);
    totalsize += size;
}


static void write_int (int fh, int *ip, int count)
{
    int size;

    size = count * sizeof (int);
    if (size != write (fh, ip, size))
        error_handler ("error writing");

    totalsize += size;
}

/******/
