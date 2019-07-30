#include "negf_prototypes.h"
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
#include "init_var.h"
#include "LCR.h"
#include "prototypes_on.h"


void is_vloc_state_overlap (STATE *states)
{

    int ion, ista;
    double r, r1, r2;
    SPECIES *sp;
    ION *iptr;

    if(vloc_state_overlap_or_not == NULL)
        my_malloc(vloc_state_overlap_or_not, pct.n_ion_center_loc * ct.num_states, char);

    for (ion = 0; ion < pct.n_ion_center_loc; ion++)
    {
        iptr = &Atoms[ion];
        sp = &ct.sp[iptr->species];
        r1 = sp->lradius;

        for (ista = 0; ista < ct.num_states; ista++)
        {

            r = minimage1 (Atoms[ion].crds, states[ista].crds);
            r2 = states[ista].radius;

            if (r < (r1 + r2))
                vloc_state_overlap_or_not[ion * ct.num_states + ista] = 1;
            else
                vloc_state_overlap_or_not[ion * ct.num_states + ista] = 0;


        } 
    }    

}
