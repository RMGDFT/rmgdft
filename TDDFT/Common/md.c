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
#include "svnrev.h"



/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
//CONTROL ct;

/* PE control structure which is also declared extern in main.h */
//PE_CONTROL pct;


int main(int argc, char **argv)
{



    time_t tt;
    char *timeptr;

    double *Hmatrix, *Smatrix;

    time(&tt);
    timeptr = ctime(&tt);
    ct.time0 = my_crtc();

    ct.images_per_node = 1;
    init_IO(argc, argv);


    my_barrier();

    
    int n2 = ct.num_states * ct.num_states;
    my_malloc( Hmatrix, n2, double );
    my_malloc( Smatrix, n2, double );

    /*  Begin to do the real calculations */
    init_TDDFT();

    for(ct.scf_steps = 0; ct.scf_steps < 5; ct.scf_steps++)
    {
        get_cholesky_real(matB);

        get_dm_diag_p(states, l_s, mat_X, Hij);

        write_eigs(states);

        update_TDDFT(mat_X);

    }

    MPI_Finalize();

    return 0;
}



