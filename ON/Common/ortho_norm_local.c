/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"

#define 	MAX_ORBIT_ON_ION 	24

void state_minus_state (STATE *state1, STATE *state2, double factor);


void ortho_norm_local (STATE *states)
{
    int ione = 1, i, j, k;
    int num_state_on_this_center;
    double norm;
    double overlap[MAX_ORBIT_ON_ION];
    double scale, x, y, z, center_dis;
    

    if (pct.gridpe == 0)
        printf ("\n LOCAL ORTHONORMALIZATION  ");


    for (i = ct.state_begin; i < ct.state_end; i++)
    {
        norm = QMD_ddot (states[i].size, states[i].psiR, ione, states[i].psiR, ione);
        norm = 1.0 / (sqrt (norm));
        dscal (&states[i].size, &norm, states[i].psiR, &ione);

        num_state_on_this_center = 0;
        for (j = i+1; j < ct.state_end; j++)
        {

            x = (states[i].crds[0] - states[j].crds[0]);
            y = (states[i].crds[1] - states[j].crds[1]);
            z = (states[i].crds[2] - states[j].crds[2]);
            center_dis = x*x + y*y + z*z;
            if (center_dis < 1.0e-6)   // two orbital has same center
            {
                num_state_on_this_center++;
                for (k = 1; k <= num_state_on_this_center; k++)
                {
                    overlap[k] = QMD_ddot (states[j].size, states[j].psiR, ione, states[j-k].psiR, ione);
                }
                for (k = 1; k <= num_state_on_this_center; k++)
                {


                    scale = -overlap[k];
                    daxpy(&states[j].size, &scale, states[j-k].psiR, &ione, states[j].psiR, &ione);
                }

                norm = QMD_ddot (states[j].size, states[j].psiR, ione, states[j].psiR, ione);
                norm = 1.0 / (sqrt (norm));
                dscal (&states[j].size, &norm, states[j].psiR, &ione);

            }
        }
        i += num_state_on_this_center;
    }

    for (j = ct.state_begin; j < ct.state_end; j++)
    {
        norm = 1.0 / (sqrt (get_vel()));
        dscal (&states[j].size, &norm, states[j].psiR, &ione);
    }

}


