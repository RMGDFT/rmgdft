/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/rmg_timings.c *****
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
 *   void write_timings(void)
 *   some other subroutines to define clock for different platform
 * INPUTS
 *   nothing
 * OUTPUT
 *   print the timeings in the standard output.
 * PARENTS
 *   too many
 * CHILDREN
 *   no
 * SOURCE
 */


#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include <time.h>
#include <stdio.h>

#include <sys/time.h>

rmg_double_t my_crtc (void)
{
    rmg_double_t val, val1;
    struct timeval t1;
    gettimeofday (&t1, NULL);
    val1 = (rmg_double_t) t1.tv_usec;
    val1 /= 1000000.0;
    val = (rmg_double_t) t1.tv_sec;
    val += val1;
    return val;
}


