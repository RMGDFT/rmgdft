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
#include "svnrev.h"
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <unordered_map>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
//#include "main.h"
#include "prototypes_on.h"
#include "transition.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"



#include "twoParts.h"
#include "pmo.h"
#include "cei.h"


#include "init_var.h"


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

std::unordered_map<std::string, InputKey *> ControlMap;

COMPASS potentialCompass, chargeDensityCompass;


void ReadBranchNEGF(char *cfile, CONTROL& lc, complex_energy_integral& cei, COMPASS& potcompass, COMPASS& rhocompass);

int main (int argc, char **argv)
{


    RmgTimer *RT = new RmgTimer("1-TOTAL");


    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;

 
#if GPU_ENABLED
    magma_init();
#endif

    ct.images_per_node = 1;
    InitIo(argc, argv, ControlMap);

    ReadBranchNEGF(ct.cfile, ct, cei, potentialCompass, chargeDensityCompass);
    allocate_states();
    ReadOrbitals (ct.cfile, states, state_to_ion, pct.img_comm);
    get_state_to_proc(states);

    my_barrier ();

    /*  Begin to do the real calculations */
    run (states, states1);


    delete(RT);

    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(Rmg_G, ct.logname, ct.scf_steps);


    MPI_Finalize ();


    return 0;
}                               /*   end main */


/******/
