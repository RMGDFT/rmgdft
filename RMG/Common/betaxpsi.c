/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void betaxpsi (STATE * states)
{

    int kpt;

    for (kpt = 0; kpt < ct.num_kpts; kpt++)
        betaxpsi1 (&states[kpt * ct.num_states], kpt);
}
