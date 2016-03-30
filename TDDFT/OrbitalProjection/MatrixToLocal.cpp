/************************** SVN Revision Information **************************
 **    $Id: tri_to_local.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"



void MatrixToLocal (STATE *states_distribute, double * A_glob, double * A_local)
{


    int j,  k;
    int jj, kk;

    for(j =0; j < pct.num_local_orbit; j++)
        for(k =0; k < pct.num_local_orbit; k++)
        {
            jj = states_distribute[j].istate;
            kk = states_distribute[k].istate;

            A_local[k * pct.num_local_orbit + j] = A_glob[kk * ct.num_states + jj] ;


        }
}


