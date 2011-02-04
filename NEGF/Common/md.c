/************************** SVN Revision Information **************************
 **    $Id: md.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/md.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   int main(int argc, char **argv)
 *   Main program
 *   Read-in all informations, structures, pseudopentials, etc. 
 *   Then enters the main driver loop. 
 * INPUTS
 *   when we run it, we need to give the input control 
 *   file name in the first argument
 *   for example, md in.diamond8
 * OUTPUT
 *   nothing
 * PARENTS
 *   This is grand-grand-....
 * CHILDREN
 *   run.c
 * SEE ALSO
 *   md.h for structure defination
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "md.h"
#include "pmo.h"
#include "svnrev.h"


STATE states[MAX_STATES];
STATE states1[MAX_STATES];

/* Main control structure which is declared extern in md.h so any module */
/* may access it.					                 */
CONTROL ct;

complex_energy_integral cei;
parallel_matrix_operation pmo;

/* Main control for two parts  --Qingzhong */
CONTROL part1;
CONTROL part2;
PE_CONTROL pct1;
PE_CONTROL pct2;

/* PE control structure which is also declared extern in md.h */
PE_CONTROL pct;

int mpi_nprocs;
int mpi_myrank;

int total_mem = 0;

int main (int argc, char **argv)
{

    char *timeptr;
    time_t tt;

    /* initialize the MPI */
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_myrank);
    pct.thispe = mpi_myrank;

    /* Read in the name of the control file from the command line */
    strcpy (ct.cfile, argv[1]);

    time (&tt);
    timeptr = ctime (&tt);

    if (pct.thispe == 0)
    {
        printf ("\n  Code Revision %d, Last change on %s", SVN_REV, SVN_REVDATE); 
        printf ("\n  Run started at %s", timeptr);
        printf ("\n  Total procs: %d \n", mpi_nprocs);

    }


    init_IO();

    /* Read in our control information */
    read_control ();

    read_trans (&cei);

    read_LCR ();

    /* Read in our pseudopotential information */
    read_pseudo ();

    my_barrier ();

    /*  Begin to do the real calculations */
    run (states, states1);

    {
        char tag[10];
        sprintf( tag, "%d", ct.runflag );
        my_alloc_report( tag );
    }

    MPI_Finalize ();


    return 0;
}                               /*   end main */


/******/
