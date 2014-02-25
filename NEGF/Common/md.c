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
 *   main.h for structure defination
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "init_var.h"
#include "twoParts.h"
#include "pmo.h"
#include "svnrev.h"


/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

complex_energy_integral cei;
parallel_matrix_operation pmo;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

int mpi_nprocs;
int mpi_myrank;

int total_mem = 0;

/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;


int main (int argc, char **argv)
{

    char *timeptr;
    time_t tt;



    time (&tt);
    timeptr = ctime (&tt);
    ct.time0 = my_crtc();


    ct.images_per_node = 1;
    init_IO(argc, argv);

    /* Read in our control information */
//    read_control ();

    read_trans (&cei);

    read_LCR ();

    /* Read in our pseudopotential information */
//    read_pseudo ();

    my_barrier ();

    /*  Begin to do the real calculations */
    run (states, states1, states_distribute);


    MPI_Finalize ();


    return 0;
}                               /*   end main */


/******/
