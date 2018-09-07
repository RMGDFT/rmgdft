/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "transition.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

extern std::vector<ORBITAL_PAIR> OrbitalPairs;


void ThetaPhiBlock(int pair_start, int pair_end, double *work_theta)
{


    for(int ipair = pair_start; ipair < pair_end; ipair++)
    {
        ORBITAL_PAIR onepair = OrbitalPairs[ipair];
        int st1 = onepair.orbital1;
        int st2 = onepair.orbital2;
        int st11 = st1-ct.state_begin;

        double temp = work_theta[st11 * ct.num_states + st2];
        ThetaPhi(st1, st2, temp, states[st2].psiR, states1[st1].psiR, 0, states, &onepair);
    }
}
