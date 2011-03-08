/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void get_inverse_block (doublecomplex *Hii, doublecomplex *Gii, int *ipiv, int nn)
{
/*  direct inverse of a small matrix  Hii  */

    int i, j, info;

    for (i = 0; i < nn * nn; i++)
    {
        Gii[i].r = 0.0;
        Gii[i].i = 0.0;
    }

    for (i = 0; i < nn; i++)
    {
        Gii[i + i * nn].r = 1.0;
    }

    ZGESV (&nn, &nn, Hii, &nn, ipiv, Gii, &nn, &info);

    if (info != 0)
    {
        printf ("get_inverse_block.c: error in ZGESV with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }


}
