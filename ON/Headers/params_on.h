/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/params.h *****
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
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */

/*

                            params.h


    Compile time parameters for md program.


*/


/*************************************************************
 *************************************************************

    All compile time parameters are now set in the top-level
    Makefile. Don't touch anything in this file unless you
    know what you are doing.

 *************************************************************
 *************************************************************/
#include "arch.h"

int PE_X, PE_Y, PE_Z;
#define  GAMMA_PT    1

/* This is the maximum number of non-local projectors for any species
 * just an upper limit on dynamically allocated data structures. */

/* Maximum cpu grids, now spin enabled -> 2.
 *  * This will need to be dynamic if we parallelize over k-points */

#define         MAX_OPTS  16

/* Maximum number of states -- just an upper limit on 
 * dynamically allocated data structures. */
#define         MAX_STATES  30000


/* Maximum number of species -- just an upper limit on
 * dynamically allocated data structures. */
#define MAX_TRADE_IMAGES 5


/* Maximum number of ions -- just an upper limit on
 * dynamically allocated data structures. */
#define         IONS_PER_PE     MAX_IONS/4


/* Number of k-points -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_KPTS           32




/* Size of the linear interpolation grid for the potentials */
#define         MAX_LOCAL_LIG       (9000)

/* Size of the linear interpolation grid for the qfunction */
/* Max. number of states localized on the same ion */
#define MAX_LOC_ST      20

/* Maximum number of images for finite difference routines */
#define MAX_TRADE_IMAGES 5


/******/
