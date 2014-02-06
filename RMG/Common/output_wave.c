/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/output_wave.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2001  Serge Nakhmanson, Wenchang Lu, Jerzy Bernholc
 * FUNCTION
 *   void output_wave(STATE *states, int kpt, int fhand)
 *   output the wavefunction for bandstructure calculations
 * INPUTS
 *   states: pointer to wavefunctions
 *   kpt:    k-point index
 *   fhand:  file name 
 * OUTPUT
 *   write to a file 
 * PARENTS
 *   bandstructure.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <unistd.h>
#include <errno.h>
#include "main.h"
#include "common_prototypes.h"

void output_wave (STATE * states, int ik, int fhand)
{
    int status, size, is;

    /* Wait until everyone gets here */
    my_barrier ();

    status = write (fhand, &ik, sizeof (int));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[0], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[1], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[2], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kweight, sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);


    size = (GAMMA_PT) ? get_P0_BASIS() : 2 * get_P0_BASIS();

    for (is = 0; is < ct.num_states; is++)
    {
        status = write (fhand, &states[is].eig, sizeof (rmg_double_t));
        if (status == -1)
            error_handler ("Unable to write wave state. ERRNO is %d.", errno);
        status = write (fhand, states[is].psiR, size * sizeof (rmg_double_t));
        if (status == -1)
            error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    }

}                               /* end write_data */

/******/
