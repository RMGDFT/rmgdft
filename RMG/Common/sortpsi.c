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



#include "common_prototypes.h"
#include "main.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>


void sortpsi (STATE * states)
{

    int state, n, incx, idx1, koffset, koffset2, loffset;
    rmg_double_t t1;
    STATE *sp, *sp1;
    rmg_double_t *tmp_psi1R, *tmp_psi2R;
    rmg_double_t *tmp_psi1I, *tmp_psi2I;
    ION *iptr;

    my_malloc (tmp_psi1R, 2 * get_P0_BASIS(), rmg_double_t);
    my_malloc (tmp_psi2R, 2 * get_P0_BASIS(), rmg_double_t);
    tmp_psi1I = NULL;
    tmp_psi2I = NULL;
#if !GAMMA_PT
    tmp_psi1I = tmp_psi1R + get_P0_BASIS();
    tmp_psi2I = tmp_psi2R + get_P0_BASIS();
#endif


    n = get_P0_BASIS();
#if !GAMMA_PT
    n *= 2;
#endif

    incx = 1;

    for (state = 0; state < ct.num_states - 1; state++)
    {

        sp = &states[state];
        sp1 = &states[state + 1];
    
	koffset = sp->kidx * ct.num_ions * ct.num_states * ct.max_nl;
	koffset2 = sp->kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl;
    

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

                t1 = sp->oldeig[0];
                sp->oldeig[0] = sp1->oldeig[0];
                sp1->oldeig[0] = t1;

                t1 = sp->occupation[0];
                sp->occupation[0] = sp1->occupation[0];
                sp1->occupation[0] = t1;

		for (idx1 = 0; idx1 < pct.num_nonloc_ions; idx1++)
                {

		    /* For localized <beta|psi>, there is offset due to both k-point and ion*/
		    loffset = koffset2 + idx1 * ct.num_states * ct.max_nl;
		    
		    my_swap(&pct.newsintR_local[loffset + state * ct.max_nl], 
			    &pct.newsintR_local[loffset + (state + 1) * ct.max_nl], ct.max_nl);
			
		    my_swap(&pct.oldsintR_local[loffset + state * ct.max_nl], 
			    &pct.oldsintR_local[loffset + (state + 1) * ct.max_nl], ct.max_nl);
                        
#if !GAMMA_PT
		    my_swap(&pct.newsintI_local[loffset + state * ct.max_nl], 
			    &pct.newsintI_local[loffset + (state + 1) * ct.max_nl], ct.max_nl);
			
		    my_swap(&pct.oldsintI_local[loffset + state * ct.max_nl], 
			    &pct.oldsintI_local[loffset + (state + 1) * ct.max_nl], ct.max_nl);

#endif
                }               /* end for  idx1 */


            }                   /* end if */

        }                       /* end if */

    }                           /* end for */


    my_free (tmp_psi2R);
    my_free (tmp_psi1R);

}                               /* end sort_psi */

/******/
