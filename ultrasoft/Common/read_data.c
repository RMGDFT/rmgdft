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

    /* wait until everybody gets here */
    my_barrier ();


    /* Make the new output file name */
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


    if (pct.gridpe == 0)
    {
        printf ("read_data: psi grid = %d %d %d\n", grid[0], grid[1], grid[2]);
        printf ("read_data: pe grid = %d %d %d\n", pe[0], pe[1], pe[2]);
        printf ("read_data: grid_size  = %d\n", grid_size);
        printf ("read_data: fine = %d %d %d\n", fine[0], fine[1], fine[2]);
        printf ("read_data: fgrid_size = %d\n", fgrid_size);
    }




    /* read wavefunction info */
    read_int (fhand, &gamma, 1);
    if (gamma != GAMMA_PT)
        error_handler ("Wrong gamma data");


    read_int (fhand, &nk, 1);
    if (nk != ct.num_kpts && ct.forceflag != BAND_STRUCTURE)    /* bandstructure calculation */
        error_handler ("Wrong number of k points");


    read_int (fhand, &ns, 1);
    if (ns != ct.num_states)
        error_handler ("Wrong number of states");

    if (pct.gridpe == 0)
    {
        printf ("read_data: gamma = %d\n", gamma);
        printf ("read_data: nk = %d\n", ct.num_kpts);
        printf ("read_data: ns = %d\n", ct.num_states);
    }



    /* read the hartree potential */
    read_double (fhand, vh, fgrid_size);
    if (pct.gridpe == 0)
        printf ("read_data: read 'vh'\n");

    /* read density */
    read_double (fhand, rho, fgrid_size);
    if (pct.gridpe == 0)
        printf ("read_data: read 'rho'\n");

    /* read Vxc */
    read_double (fhand, vxc, fgrid_size);
    if (pct.gridpe == 0)
        printf ("read_data: read 'vxc'\n");






    /* read state occupations */
    {
        STATE *sp;
        REAL *occ;
        my_malloc (occ, nk * ns, REAL);
        read_double (fhand, occ, nk * ns);

        if (pct.gridpe == 0)
            printf ("read_data: read 'occupations'\n");

        if (ct.override_occ != 1 && ct.forceflag != BAND_STRUCTURE)
        {
            REAL occ_total = 0.0;
            sp = states;
            for (ik = 0; ik < nk; ik++)
                for (is = 0; is < ns; is++)
                {

                    occ_total += (sp->occupation = occ[ik * ns + is]);
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
                        sp->occupation *= fac;
                        sp++;
                    }
            }

        }
        my_free (occ);

    }




    /* read state eigenvalues, not needed really */
    {


        STATE *sp;

        sp = states;
        for (ik = 0; ik < nk; ik++)
            for (is = 0; is < ns; is++)
            {

                read_double (fhand, &sp->eig, 1);
                sp++;
            }

        if (pct.gridpe == 0)
            printf ("read_data: read 'eigenvalues'\n");

    }


    /* read wavefunctions */
    {
        int wvfn_size = (gamma) ? grid_size : 2 * grid_size;
        STATE *sp;

#ifdef SMP
        /* get temporary buffers */
        REAL *work1, *work2;
        my_malloc (work1, wvfn_size, REAL);
        work2 = (gamma) ? NULL : work1 + grid_size;
#endif

        sp = states;
        for (ik = 0; ik < nk; ik++)
        {
            for (is = 0; is < ns; is++)
            {

#ifdef SMP
                read_double (fhand, work1, wvfn_size);
                scatter_psi (work1, work2, sp, 0);
#else
                read_double (fhand, sp->psiR, wvfn_size);
#endif

                sp++;
            }
            /*  for calculating band structures, 
               only read-in one wave function (?) */
            if (ct.forceflag == BAND_STRUCTURE)
                return;
        }

#ifdef SMP
        /* Release memory */
        my_free (work1);
#endif
        if (pct.gridpe == 0)
            printf ("read_data: read 'wfns'\n");

    }






    /* read number of ions */
    read_int (fhand, &na, 1);
    if (na != ct.num_ions)
        error_handler ("Wrong number of ions");


    /* read current ionic cartesian positions */
    {
        REAL r[3];
        for (ia = 0; ia < na; ia++)
        {
            read_double (fhand, r, 3);

            if (ct.override_current)
                for (i = 0; i < 3; i++)
                    ct.ions[ia].crds[i] = r[i];

        }


        /* read current ionic crystal positions */
        for (ia = 0; ia < na; ia++)
        {
            read_double (fhand, r, 3);

            if (ct.override_current)
            {
                for (i = 0; i < 3; i++)
                    ct.ions[ia].xtal[i] =
                        (r[i] < 0) ? (r[i] + 1.0) : ((r[i] > 1.0) ? (r[i] - 1.0) : r[i]);

                to_cartesian (ct.ions[ia].xtal, ct.ions[ia].crds);
            }

        }

        /* Overwrite the initial positions with current positions, 
         * may be useful for constraint dynamics  */

        /* read original ionic cartesian positions */
        for (ia = 0; ia < na; ia++)
        {
            read_double (fhand, &ct.ions[ia].icrds[0], 3);

            if (ct.override_initial)
                for (i = 0; i < 3; i++)
                    ct.ions[ia].icrds[i] = ct.ions[ia].crds[i];
        }


        /* read original ionic crystal positions */
        for (ia = 0; ia < na; ia++)
        {
            read_double (fhand, &ct.ions[ia].ixtal[0], 3);

            if (ct.override_initial)
                for (i = 0; i < 3; i++)
                    ct.ions[ia].ixtal[i] = ct.ions[ia].xtal[i];
        }


        /* read ionic velocities */
        for (ia = 0; ia < na; ia++)
        {
            read_double (fhand, &ct.ions[ia].velocity[0], 3);
        }

        /* read forces pointer */
        read_int (fhand, &ct.fpt[0], 4);

        /* read ionic forces */
        for (ia = 0; ia < na; ia++)
            read_double (fhand, &ct.ions[ia].force[0][0], 3 * 4);


        /* read Nose positions,velocities, masses and forces from the file */
        read_double (fhand, &ct.nose.xx[0], 10);
        read_double (fhand, &ct.nose.xv[0], 10);
        read_double (fhand, &ct.nose.xq[0], 10);
        read_double (fhand, &ct.nose.xf[0][0], 4 * 10);

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
