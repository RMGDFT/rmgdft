/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"


/* Reads and parses the input control file */
void change_states_crds(STATE * states)
{

    int ion, ist;

    /* Coordinates and species type for each ion. */
    for (ist = 0; ist < ct.num_states; ist++)
    {

        ion = state_to_ion[ist];

        if (ct.ions[ion].movable)
        {

            states[ist].crds[0] = ct.ions[ion].crds[0];
            states[ist].crds[1] = ct.ions[ion].crds[1];
            states[ist].crds[2] = ct.ions[ion].crds[2];


        }                       /* end if */

    }


}
