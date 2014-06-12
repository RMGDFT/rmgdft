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
#include "prototypes_on.h"
#include "init_var.h"
#include "svnrev.h"



/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

int mpi_nprocs;
int mpi_myrank;


/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;


int main(int argc, char **argv)
{


    void *RT = BeginRmgTimer("1-TOTAL");
    ct.images_per_node = 1;
    init_IO(argc, argv);



    init_states();
    my_barrier();

    /*  Begin to do the real calculations */
    run(states, states1);


    EndRmgTimer(RT);


    if(pct.imgpe == 0) fclose(ct.logfile);
    CompatRmgTimerPrint(ct.logname, ct.scf_steps);

    MPI_Finalize();

    return 0;
}


