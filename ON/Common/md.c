/************************** SVN Revision Information **************************
 **    $Id$    **
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
#include "svnrev.h"



/* Main control structure which is declared extern in md.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in md.h */
PE_CONTROL pct;

int mpi_nprocs;
int mpi_myrank;


int main(int argc, char **argv)
{



    time_t tt;
    char *timeptr;


    time(&tt);
    timeptr = ctime(&tt);
    ct.time0 = my_crtc();

    ct.images_per_node = 1;
    init_IO(argc, argv);


    /* Read in the name of the control file from the command line */
//    strcpy(ct.cfile, argv[1]);

    /* Read in our control information */
 //   read_control();

    /* Read in our pseudopotential information */
//    read_pseudo();

    my_barrier();

    /*  Begin to do the real calculations */
    run(states, states1);

    //my_alloc_report( "order-n" );

    write_timings();
    MPI_Finalize();

    return 0;
}


