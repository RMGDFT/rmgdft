/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                        get_normKB.c


*/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"



/* Set up Kleinman-Bylander normalization coefficients */

void get_normKB(SPECIES * sp, double *pd)
{
    int pcount, ip;


    pcount = 0;
    for (ip = 0; ip < sp->num_potentials; ip++)
    {

        if (sp->lval[ip] != sp->local)
        {

            switch (sp->lval[ip])
            {

            case S_STATE:
                pd[pcount] = sp->kbnorm[ip];
#if DEBUG
                if (pct.thispe == 0)
                    printf(" kbnorm[%d]= %f\n", ip, sp->kbnorm[ip]);
#endif
                pcount++;
                break;

            case P_STATE:
                pd[pcount] = 3.0 * sp->kbnorm[ip];
                pd[pcount + 1] = 3.0 * sp->kbnorm[ip];
                pd[pcount + 2] = 3.0 * sp->kbnorm[ip];
                pcount += 3;
                break;

            case D_STATE:
                pd[pcount] = 5.0 * sp->kbnorm[ip];
                pd[pcount + 1] = 5.0 * sp->kbnorm[ip];
                pd[pcount + 2] = 5.0 * sp->kbnorm[ip];
                pd[pcount + 3] = 5.0 * sp->kbnorm[ip];
                pd[pcount + 4] = 5.0 * sp->kbnorm[ip];
                pcount += 5;
                break;
            default:
                error_handler("Angular momentum state not programmed");


            }                   /* end switch */

        }                       /* end if */

    }                           /* end for */

}
