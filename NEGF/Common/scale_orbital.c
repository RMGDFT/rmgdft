/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

    read orbitals which has overlap with this processor domain and map it on
*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void scale_orbital (STATE * states, STATE *states_distribute)
{


    int st2, idx;
    for (st2 = 0; st2 < pct.num_local_orbit; st2++)
    {
        for(idx = 0; idx < get_P0_BASIS(); idx++)
            states_distribute[st2].psiR[idx] *= 1.0/sqrt(get_vel());
    }   

    for(st2 = ct.state_begin; st2 < ct.state_end; st2++)
    for(idx = 0; idx < states[st2].size; idx++)
        states[st2].psiR[idx] *= 1.0/sqrt(get_vel());
}   


