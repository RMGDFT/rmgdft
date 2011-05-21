/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/scatter_psi.c *****
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
 *   void scatter_psi(REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid)
 *   For smp mode, this function is used to scatter an orbital from a local
 *   memory array to remote memory
 *   For mpi mode, copy tmp_psiR to sp->psiR and tmp_psiI to sp->psiI 
 * INPUTS
 *   tmp_psiR: real part of wave function
 *   tmp_psiI: imaginary part of wave function
 *   tid: Thread ID
 * OUTPUT
 *   sp:  points to orbital structure STATE (see main.h)
 * PARENTS
 *   init_wf.c init_wflcao.c mg_eig_state.c norm_psi.c sortpsi.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <stdlib.h>
#include "main.h"



void scatter_psi (REAL * tmp_psiR, REAL * tmp_psiI, STATE * sp, int tid)
{

    int idx;
    if (tmp_psiR != NULL)
    {
        for (idx = 0; idx < sp->pbasis; idx++)
            sp->psiR[idx] = tmp_psiR[idx];
    }                           /* end if */

    if (tmp_psiI != NULL)
    {
        for (idx = 0; idx < sp->pbasis; idx++)
            sp->psiI[idx] = tmp_psiI[idx];
    }                           /* end if */

}                               /* end scatter_psi */


/******/
