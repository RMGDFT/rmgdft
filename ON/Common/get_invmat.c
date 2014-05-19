/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

/* 
   	output: matB = (matB)^-1
*/



#include "my_scalapack.h"


void get_invmat(double *matB)
{
    int numst = ct.num_states;
    int info;
    double time1, time2;
    char uplo = 'l';

    int ione = 1;
    int myrow;
    _fcd char_fcd1;

    time1 = my_crtc();


    myrow = pct.myrow;

    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {
        /* Compute the Cholesky decomposition of matB */
        char_fcd1 = &uplo;
        PSPOTRF(char_fcd1, &numst, matB, &ione, &ione, pct.desca, &info);
        if (info != 0)
        {
            printf(" PSPOTRF in get_invmat.c, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }
        /* now get the inverse */
        PSPOTRI(char_fcd1, &numst, matB, &ione, &ione, pct.desca, &info);
        if (info != 0)
        {
            printf(" PSPOTRI, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }
    }

    my_barrier();


    time2 = my_crtc();
    rmg_timings(INVMATB_TIME, (time2 - time1));
}
