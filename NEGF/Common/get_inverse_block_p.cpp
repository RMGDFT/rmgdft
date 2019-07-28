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

void get_inverse_block_p (std::complex<double> *Hii, std::complex<double> *Gii, int *ipiv, int *desca )
{
/*  direct inverse of a small matrix  Hii  */

    int i, j, info;
    int ione =1;
    int nn = desca[2];


    pmo_unitary_matrix(Gii, desca);
    pzgesv (&nn, &nn, Hii, &ione, &ione, desca, ipiv, Gii, &ione, &ione, desca, &info);

    if (info != 0)
    {
        printf ("get_inverse_block_p.c: error in pzgesv with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }


}
