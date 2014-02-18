/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main_on.h"


void normalize_orbits(STATE * states)
{
    int st, i;
    REAL sum;
    REAL tmp;
    int ione = 1;

    if (pct.gridpe == 0)
        printf("\n STATE  NORMALIZATION  ");

    for (st = ct.state_begin; st < ct.state_end; st++)
    {

        sum = 0.0;
        for (i = 0; i < states[st].size; i++)
            sum += states[st].psiR[i] * states[st].psiR[i];

        sum = sqrt(sum);

        for (i = 0; i < states[st].size; i++)
            states[st].psiR[i] /= sum;
        tmp = 1. / sqrt(ct.vel);
        sscal(&states[st].size, &tmp, states[st].psiR, &ione);

    }

}
