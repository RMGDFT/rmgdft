/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/sortpsi.c *****
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
 *   void sortpsi(STATE *states)
 *   Sorts the input orbitals by eigenvalue (ascending)
 *   It presumes that they are already nearly sorted and merely does
 *   nearest neighbor swapping.
 *   It never swaps occupied and unoccupied states.
 *
 *   This routine only sorts the states associated with one
 *   k-point at a time 
 * INPUTS
 *   states: point to orbital structure (see main.h)
 * OUTPUT
 *   states are updated
 * PARENTS
 *   moldyn.c scf.c
 * CHILDREN
 *   gather_psi.c scatter_psi.c
 * SOURCE
 */



#include "main.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>


void sortpsi (STATE * states)
{

    int state, n, incx, idx1, idx2;
    REAL t1;
    STATE *sp, *sp1;
    REAL *tmp_psi1R, *tmp_psi2R;
    REAL *tmp_psi1I, *tmp_psi2I;
    ION *iptr;

    my_malloc (tmp_psi1R, 2 * P0_BASIS, REAL);
    my_malloc (tmp_psi2R, 2 * P0_BASIS, REAL);
    tmp_psi1I = NULL;
    tmp_psi2I = NULL;
#if !GAMMA_PT
    tmp_psi1I = tmp_psi1R + P0_BASIS;
    tmp_psi2I = tmp_psi2R + P0_BASIS;
#endif


    n = P0_BASIS;
#if !GAMMA_PT
    n *= 2;
#endif

    incx = 1;

    for (state = 0; state < ct.num_states - 1; state++)
    {

        sp = &states[state];
        sp1 = &states[state + 1];

        if (sp->eig[0] > sp1->eig[0])
        {

            if (((sp->occupation[0] > 0.1) && (sp1->occupation[0] > 0.1))
                || ((sp->occupation[0] < 0.1) && (sp1->occupation[0] < 0.1)))
            {
                gather_psi (tmp_psi1R, tmp_psi1I, sp, 0);
                gather_psi (tmp_psi2R, tmp_psi2I, sp1, 0);
                scatter_psi (tmp_psi1R, tmp_psi1I, sp1, 0);
                scatter_psi (tmp_psi2R, tmp_psi2I, sp, 0);

                t1 = sp->eig[0];
                sp->eig[0] = sp1->eig[0];
                sp1->eig[0] = t1;

                t1 = sp->occupation[0];
                sp->occupation[0] = sp1->occupation[0];
                sp1->occupation[0] = t1;

                for (idx1 = 0; idx1 < ct.num_ions; idx1++)
                {
                    iptr = &ct.ions[idx1];

                    for (idx2 = 0; idx2 < ct.max_nl; idx2++)
                    {

                        t1 = iptr->newsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            state * ct.max_nl + idx2];
                        iptr->newsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       state * ct.max_nl + idx2] =
                            iptr->newsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                           (state + 1) * ct.max_nl + idx2];
                        iptr->newsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       (state + 1) * ct.max_nl + idx2] = t1;

                        t1 = iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            state * ct.max_nl + idx2];
                        iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       state * ct.max_nl + idx2] =
                            iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                           (state + 1) * ct.max_nl + idx2];
                        iptr->oldsintR[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       (state + 1) * ct.max_nl + idx2] = t1;

#if !GAMMA_PT
                        t1 = iptr->newsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            state * ct.max_nl + idx2];
                        iptr->newsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       state * ct.max_nl + idx2] =
                            iptr->newsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                           (state + 1) * ct.max_nl + idx2];
                        iptr->newsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       (state + 1) * ct.max_nl + idx2] = t1;

                        t1 = iptr->oldsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            state * ct.max_nl + idx2];
                        iptr->oldsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       state * ct.max_nl + idx2] =
                            iptr->oldsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                           (state + 1) * ct.max_nl + idx2];
                        iptr->oldsintI[sp->kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                       (state + 1) * ct.max_nl + idx2] = t1;

#endif
                    }           /* end for idx2 */
                }               /* end for  idx1 */


            }                   /* end if */

        }                       /* end if */

    }                           /* end for */


    my_free (tmp_psi2R);
    my_free (tmp_psi1R);

}                               /* end sort_psi */







/* This doesn't work yet */
#if SMP

void sort_psi_smp (SCF_THREAD_CONTROL * s)
{

    int idx, state, n, incx;
    REAL t1;
    REAL *rptr1, *rptr2;

    n = s->numpt;
    incx = 1;

    for (state = 0; state < ct.num_states - 1; state++)
    {

        if (s->eigs[state] > s->eigs[state + 1])
        {

            if (((s->occs[state] > 0.1) && (s->occs[state + 1] > 0.1))
                || ((s->occs[state] < 0.1) && (s->occs[state + 1] < 0.1)))
            {

                t1 = s->eigs[state];
                s->eigs[state] = s->eigs[state + 1];
                s->eigs[state + 1] = t1;

                t1 = s->occs[state];
                s->occs[state] = s->occs[state + 1];
                s->occs[state + 1] = t1;

                rptr1 = &s->base_mem[state * s->lda];
                rptr2 = &s->base_mem[(state + 1) * s->lda];
                sswap (&n, rptr1, &incx, rptr2, &incx);

            }                   /* end if */

        }                       /* end if */

    }                           /* end for */

}                               /* end sort_psi_smp */

#endif
/******/
