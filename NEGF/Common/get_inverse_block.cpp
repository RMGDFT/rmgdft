#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"

void get_inverse_block (std::complex<double> *Hii, std::complex<double> *Gii, int *ipiv, int nn)
{
/*  direct inverse of a small matrix  Hii  */

    int i,info;

    for (i = 0; i < nn * nn; i++)
    {
        Gii[i] = 0.0;
    }

    for (i = 0; i < nn; i++)
    {
        Gii[i + i * nn] = 1.0;
    }

    zgesv (&nn, &nn, (double *)Hii, &nn, ipiv, (double *)Gii, &nn, &info);

    if (info != 0)
    {
        printf ("get_inverse_block.c: error in zgesv with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }


}
