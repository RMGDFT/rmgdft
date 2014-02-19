/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*** QMD-MGDFT/main.h *****
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

                        main.h

    Header file for molecular dynamics code.



*/


/* Typedefs */
#include    "main_share.h"


#include "lbfgs.h"


/* Prototypes for function calls*/
#include "prototypes.h" 
//#include "common_prototypes.h" 

/* Some definitions needed for using compile time values for global grid sizes */
#include "fixed_dims.h"

/* Include the  library of exchange and correlation (namely, Libxc) */
/* #include "../../Common/libxc/include/xc.h"   */


/******/
