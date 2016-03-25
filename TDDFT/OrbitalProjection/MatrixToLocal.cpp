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

    for(j =0; j < ct.num_states; j++)
    {
        for(k=0; k < ct.num_states; k++)
        {
            jj = states_distribute[j].local_index;
            kk = states_distribute[k].local_index;

            if(jj >= 0 && kk >=0)
                A_local[kk * pct.num_local_orbit + jj] = A_glob[k * ct.num_states + j] ;


        }
    }

}

