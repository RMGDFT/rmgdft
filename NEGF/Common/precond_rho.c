/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "md.h"



void precond_rho (double *res)
{
    int cycles, ione = 1, nits, sbasis;
    double diag, t1, time1, time2, d1, one = 1.;
    double *sg_rho, *rho_res, *work1, *work2;
    int idx;

    error_handler("\n not programed for precond_rho ");

}
