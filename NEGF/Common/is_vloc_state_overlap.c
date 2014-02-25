/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

is_atom_overlap.c

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"


void is_vloc_state_overlap (STATE *states)
{

    int ion, ista;
    rmg_double_t r, r1, r2;
    SPECIES *sp;
    ION *iptr;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
	iptr = &ct.ions[ion];
	sp = &ct.sp[iptr->species];
	r1 = sp->lradius;

        for (ista = 0; ista < ct.num_states; ista++)
        {

            r = minimage1 (ct.ions[ion].crds, states[ista].crds);
            r2 = states[ista].radius;

            if (r < (r1 + r2))
                vloc_state_overlap_or_not[ion * ct.num_states + ista] = 1;
            else
                vloc_state_overlap_or_not[ion * ct.num_states + ista] = 0;


        } 
    }    

}
