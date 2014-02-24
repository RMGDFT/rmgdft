/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/rand0.c *****
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
 *   rmg_double_t rand0(long *idum)
 *   Portable random number generator for producing repeatable random
 *   starts across computer architectures and processor topologies.
 *   The basic random algorithm used here is based on the Minimal Standard 
 *   Generator proposed by Lewis, Goodman, and Miller in 1969 as implemented in
 *   "Numerical Recipes in C. 2nd. Edition"
 * INPUTS
 *   idum: random number
 * OUTPUT
 *   a random number is returned
 * PARENTS
 *   init_wf.c init_wflcao.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include "main.h"
#include "prototypes_on.h"


#define   IA  16807
#define   IM  2147483647
#define   AM  (1.0/IM)
#define   IQ  127773
#define   IR  2836
#define   MASK 123459876



rmg_double_t rand0(long *idum)
{


    long k;
    rmg_double_t ans;

    *idum ^= MASK;
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
        *idum += IM;
    ans = AM * (rmg_double_t) (*idum);
    *idum ^= MASK;
    return ans;


}                               /* end rand0 */

/******/
