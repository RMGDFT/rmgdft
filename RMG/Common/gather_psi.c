/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/gather_psi.c *****
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
 *   void gather_psi(rmg_double_t *tmp_psiR, rmg_double_t *tmp_psiI, STATE *sp, int tid)
 *   For smp mode, this function is used to gather an orbital that 
 *   is distributed across remote memory into a local memory array.
 *   For mpi mode, copy sp->psiR to tmp_psiR and sp->psiI to tmp_psiI
 * INPUTS
 *   sp:  points to orbital structure STATE (see main.h)
 *   tid: Thread ID
 * OUTPUT
 *   tmp_psiR: real part of wave function
 *   tmp_psiI: imaginary part of wave function
 * PARENTS
 *   get_ke.c get_milliken.c get_rho.c mg_eig_state.c nlforce_s.c
 *   nlforce_p.c nlforce_d.c norm_psi.c sortpsi.c subdiag_mpi.c
 *   subdiag_smp.c  write_data.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */

#include <stdlib.h>
#include "main.h"



void gather_psi (rmg_double_t * tmp_psiR, rmg_double_t * tmp_psiI, STATE * sp, int tid)
{

    int idx;
#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    if (tmp_psiR != NULL)
    {
        for (idx = 0; idx < sp->pbasis; idx++)
            tmp_psiR[idx] = sp->psiR[idx];
    }                           /* end if */

    if (tmp_psiI != NULL)
    {
        for (idx = 0; idx < sp->pbasis; idx++)
            tmp_psiI[idx] = sp->psiI[idx];
    }                           /* end if */

#if MD_TIMERS
    time2 = my_crtc ();
#endif
}                               /* end gather_psi */


/******/
