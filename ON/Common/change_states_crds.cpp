/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


/* Reads and parses the input control file */
void change_states_crds(STATE * states)
{

    int ion, ist;

    /* Coordinates and species type for each ion. */
    for (ist = 0; ist < ct.num_states; ist++)
    {

        ion = states[ist].atom_index;

        if (Atoms[ion].movable)
        {

            states[ist].crds[0] = Atoms[ion].crds[0];
            states[ist].crds[1] = Atoms[ion].crds[1];
            states[ist].crds[2] = Atoms[ion].crds[2];


        }                       /* end if */

    }


}
