/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/main.c *****
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
 *   int main(int argc, char **argv)
 *   Main program
 *   Perform any initializations that are required and then
 *   enters the main driver loop. It also handles checkpointing and the
 *   output of intermediate results.
 * INPUTS
 *   when we run it, we need to give the input control file name in the first argument
 *   for example, md in.diamond8
 * OUTPUT
 *   Print standard output.
 * PARENTS
 *   This is grand-grand-....
 * CHILDREN
 *   init_pe.c write_header.c write_occ.c quench.c fastrlx.c cdfastrlx.c moldyn.c
 *   dendx.c psidx.c write_data.c write_avgv.c write_avgd.c write_zstates.c get_milliken.c
 * SEE ALSO
 *   main.h for structure definition
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "main.h"

void initialize (int argc, char **argv);

void run (void);

void report (void);

void finish (void);



/* "VariableName_buff" used to receive and hold the same variable from the oppisite spin*/

/* State storage pointer (memory allocated dynamically in rd_cont */
STATE *states;  



/* Electronic Charge density */
REAL rho[FP0_BASIS], rho_buff[FP0_BASIS];   /* buff used to hold opposite spin density*/


/* Core Charge density */
REAL rhocore[FP0_BASIS];


/* Compensating charge density */
REAL rhoc[FP0_BASIS];


/* Hartree potential */
REAL vh[FP0_BASIS];

/* Nuclear local potential */
REAL vnuc[FP0_BASIS];

/* Exchange-correlation potential */
REAL vxc[FP0_BASIS];


/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

/*Other global variables will be defined here, so that they can be defined
 * in header files as extern as it should be. Otherwise compiler
 * on Origin complains that the variables are multiply defined*/


/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;

/*Compile time checks*/
#if (PE_X*PE_Y*PE_Z != NPES)
# error "Error: PE_X*PY_Y*PE_Z does not match NPES"
#endif


#if (NX_GRID % PE_X != 0)
# error "Number of grid points in X direction (NX_GRID) is not divisable by number of processors in that direction (PE_X)"
#endif
#if (NY_GRID % PE_Y != 0)
# error "Number of grid points in Y direction (NY_GRID) is not divisable by number of processors in that direction (PE_Y)"
#endif
#if (NZ_GRID % PE_Z != 0)
# error "Number of grid points in Z direction (NZ_GRID) is not divisable by number of processors in that direction (PE_Z)"
#endif

#if (FNX_GRID % PE_X != 0)
# error "Number of grid points in X direction (FNX_GRID) is not divisable by number of processors in that direction (PE_X)"
#endif
#if (FNY_GRID % PE_Y != 0)
# error "Number of grid points in Y direction (FNY_GRID) is not divisable by number of processors in that direction (PE_Y)"
#endif
#if (FNZ_GRID % PE_Z != 0)
# error "Number of grid points in Z direction (FNZ_GRID) is not divisable by number of processors in that direction (PE_Z)"
#endif


int main (int argc, char **argv)
{

	initialize (argc, argv);

        run ();

	report ();

	finish ();

    return 0;
}


void initialize(int argc, char **argv) {

    /* start the benchmark clock */
    ct.time0 = my_crtc ();
    
    /* Define a default output stream, gets redefined to a file in init_pe */
    ct.logfile = stdout;

    /* Initialize MPI, we need it for error_handler, amongst others */
    MPI_Init (&argc, &argv);
    

    /* Initialize all I/O including MPI group comms */
    /* Also reads control and pseudopotential files*/
    init_IO (argc, argv);


    /* initialize states */
    if (pct.spin_flag)
	    states=init_states_spin ();
    else
            states = init_states ();


    my_barrier ();

    /* Record the rime it took from the start of run until we hit init */
    rmg_timings ( PREINIT_TIME, my_crtc () - ct.time0, 0);

    /* Perform any necessary initializations */
    
    if (pct.spin_flag)
    	init_spin (vh, rho, rho_buff, rhocore, rhoc, states, vnuc, vxc);
    else 
    	init (vh, rho, rhocore, rhoc, states, vnuc, vxc);


  

    if (pct.imgpe == 0 )
    {

        /* Write header to stdout */
        write_header ();


        /* Write state occupations to stdout */
        write_occ (states);

    }


    /* Flush the results immediately */

    fflush (NULL);
    

    /* Wait until everybody gets here */
    MPI_Barrier(MPI_COMM_WORLD);

}

void run (void)
{

    REAL time2;

    /* Dispatch to the correct driver routine */
    switch (ct.forceflag)
    {

    case MD_QUENCH:            /* Quench the electrons */
        ct.max_rlx_steps = 0;
	if (pct.spin_flag)
	{
		fastrlx_spin (states, vxc, vh, vnuc, rho, rho_buff, rhocore, rhoc);
	}
	else
        	fastrlx (states, vxc, vh, vnuc, rho, rhocore, rhoc);
        break;

    case MD_FASTRLX:           /* Fast relax */
	if (pct.spin_flag)
	{
		fastrlx_spin (states, vxc, vh, vnuc, rho, rho_buff, rhocore, rhoc);
	}
	else
        	fastrlx (states, vxc, vh, vnuc, rho, rhocore, rhoc);
        break;

    case NEB_RELAX:           /* nudged elastic band relax */
        neb_relax (states, vxc, vh, vnuc, rho, rhocore, rhoc);
        break;

    case MD_CVE:               /* molecular dynamics */
    case MD_CVT:
    case MD_CPT:
        quench (states, vxc, vh, vnuc, rho, rhocore, rhoc);
        moldyn (states, vxc, vh, vnuc, rho, rhoc, rhocore);
        break;

    case BAND_STRUCTURE:
        bandstructure (states, vxc, vh, vnuc);
        break;
    default:
        error_handler ("Undefined MD method");

    }

    time2 = my_crtc ();

    rmg_timings (FINISH_TIME, (my_crtc () - time2), 0);

}                               /* end run */

void report ()
{

    /* write planar averages of quantities */
    if (ct.zaverage == 1)
    {
        /* output the average potential */
        write_avgv (vh, vnuc);
        write_avgd (rho);
    }
    else if (ct.zaverage == 2)
    {
        write_zstates (states);
    }


    /* If milliken population info is requested then compute and output it */
    if (ct.domilliken)
        mulliken (states);


    /*Destroy wisdom that may have been allocated previously */
    destroy_fftw_wisdom ();


    if (ct.write_memory_report)
    {
        /*Release memory for projectors first */
        finish_release_mem (states);

        /* write allocation report for each processor */
        my_alloc_report ("real-space");

    }

    /* Write timing information */
    if (pct.imgpe == 0)
    {
        write_timings ();
    }
}                               /* end report */


void finish ()
{

    MPI_Finalize ();

}                               /* end finish */


/******/
